"""Independent ANDES cross-check of the ieee39 line trip studied in src/line_trip_continuation.jl.

Run inside the andes conda environment:

    conda activate andes
    python scripts/andes_line_trip_check.py

The Barq side of the comparison comes from results_topology/line14_reinit.json, which
src/line_trip_continuation.jl writes.

The ANDES case is rebuilt from the workbook so that it contains only what Barq's loader
actually reads (Bus, PQ, PV, Slack, Line, GENCLS), because the workbook also carries
governors, exciters, stabilizers, shunts and a bus fault that Barq silently ignores. The
machine sheet holds GENROU-format columns but Barq only consumes M, ra and xd1, so it is
loaded here as GENCLS to keep both sides on the classical model.
"""

import json
import os

import numpy as np
import pandas as pd

import andes

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SOURCE_CASE = os.path.join(REPO, "cases/Fault_Cases/ieee39_fault.xlsx")
BARQ_RESULT = os.path.join(REPO, "results_topology/line14_reinit.json")
ANDES_CASE = "/tmp/andes_line_trip.xlsx"

TRIP_LINE = "Line_14"          # buses 8-9, the 14th row of the Line sheet
LOAD_SCALE = 1.2               # matches models.load.p[:] .*= 1.2
T_TRIP = 1.0
T_POST = 3.0                   # long enough to see whether the post-trip state holds
TSTEP = 5e-4                   # matches the Barq integrator step

# Barq's balance! splits each load evenly between constant P, I and Z
ZIP = 1.0 / 3.0

# Sheets Barq's loader consumes, plus Toggler for the trip. Everything else is dropped so
# that neither side has devices the other lacks.
KEEP_SHEETS = ["Bus", "PQ", "PV", "Slack", "Line", "GENCLS", "Area"]

# Columns belonging to GENROU that ANDES's GENCLS does not use
GENROU_ONLY = ["xd", "xq", "xd2", "xq1", "xq2", "Td10", "Td20", "Tq10", "Tq20"]


def build_case(include_shunt=False, machine_base="andes"):
    """Write an ANDES workbook holding only the devices Barq models, tripping TRIP_LINE.

    machine_base="andes" keeps the sheet as authored, so ANDES converts M, xd1, ra and D
    from each machine's Sn base to the system base. machine_base="barq" instead reproduces
    what _build_gencls does: the raw sheet numbers are used directly as system-base values,
    armature resistance is dropped and damping is pinned to the hardcoded d = 1.0.
    """
    book = pd.read_excel(SOURCE_CASE, sheet_name=None)

    # the machine sheet may be named GENROU (its columns are GENROU's) but Barq reads it as
    # a classical machine, so present it to ANDES as GENCLS
    machines = book.get("GENCLS", book.get("GENROU"))
    if machines is None:
        raise RuntimeError("no GENCLS or GENROU sheet in %s" % SOURCE_CASE)

    sheets = {name: book[name].copy() for name in KEEP_SHEETS if name in book}
    sheets["GENCLS"] = machines.drop(columns=GENROU_ONLY, errors="ignore").copy()
    if include_shunt and "Shunt" in book:
        sheets["Shunt"] = book["Shunt"].copy()

    if machine_base == "barq":
        gen = sheets["GENCLS"]
        bus_vn = book["Bus"].set_index("idx")["Vn"]
        # Sn = Vn = the system/bus rating makes every ANDES base conversion the identity,
        # so xd1 and M land in the equations exactly as Barq reads them off the sheet
        gen["Sn"] = 100.0
        gen["Vn"] = [float(bus_vn.loc[b]) for b in gen["bus"]]
        gen["ra"] = 0.0
        gen["xl"] = 0.0
        gen["D"] = 1.0
    elif machine_base != "andes":
        raise ValueError("machine_base must be 'andes' or 'barq'")

    sheets["PQ"]["p0"] = sheets["PQ"]["p0"] * LOAD_SCALE

    sheets["Toggler"] = pd.DataFrame([{
        "uid": 0, "idx": "Toggler_1", "u": 1, "name": "Toggler_1",
        "model": "Line", "dev": TRIP_LINE, "t": T_TRIP,
    }])

    with pd.ExcelWriter(ANDES_CASE) as writer:
        for name, df in sheets.items():
            df.to_excel(writer, sheet_name=name, index=False)
    return ANDES_CASE


def run_andes(case):
    """Solve power flow and a fixed-step TDS across the trip.

    Returns (system, power-flow voltages). The voltages must be copied before the TDS runs,
    because the integrator overwrites Bus.v.v in place.
    """
    ss = andes.load(case, setup=False, no_output=True, default_config=True)

    # constant PQ in power flow, matching Barq, then an even ZIP split for the TDS
    ss.PQ.config.pq2z = 0
    ss.PQ.config.p2p, ss.PQ.config.p2i, ss.PQ.config.p2z = ZIP, ZIP, ZIP
    ss.PQ.config.q2q, ss.PQ.config.q2i, ss.PQ.config.q2z = ZIP, ZIP, ZIP

    ss.setup()
    ss.PFlow.run()
    if not ss.PFlow.converged:
        raise RuntimeError("ANDES power flow did not converge")
    v_pf = np.array(ss.Bus.v.v).copy()

    ss.TDS.config.tf = T_TRIP + T_POST
    ss.TDS.config.tstep = TSTEP
    ss.TDS.config.fixt = 0              # let the step shrink, so a stall is the solver's verdict
    ss.TDS.config.criteria = 0          # do not abort on the stability criteria
    ss.TDS.init()
    ss.TDS.run()
    return ss, v_pf


def voltage_at(ss, t_target):
    """Bus voltage magnitudes at the stored time point nearest to and at or after t_target."""
    t = np.array(ss.dae.ts.t)
    y = np.array(ss.dae.ts.y)
    idx = int(np.argmax(t >= t_target)) if (t >= t_target).any() else len(t) - 1
    return t[idx], y[idx, ss.Bus.v.a]


def compare(barq, machine_base):
    """Run one ANDES configuration and report it against the Barq tracked solution."""
    buses = barq["buses"]
    ss, v_pf = run_andes(build_case(include_shunt=False, machine_base=machine_base))
    bus_ids = list(ss.Bus.idx.v)
    col = {b: i for i, b in enumerate(bus_ids)}

    t_pre, v_pre = voltage_at(ss, T_TRIP - TSTEP)    # last point before the trip
    t_post, v_post = voltage_at(ss, T_TRIP + TSTEP)  # first point after the trip

    print("\n" + "=" * 66)
    print("machine parameters: %s"
          % ("raw sheet values, as Barq reads them" if machine_base == "barq"
             else "converted from each machine's Sn to the system base"))
    print("=" * 66)
    print("  xd1 = %s" % np.round(ss.GENCLS.xd1.v, 4))
    print("  M   = %s" % np.round(ss.GENCLS.M.v, 2))
    print("  TDS reached t = %.4f of %.4f (pre sample %.4f, post sample %.4f)"
          % (ss.dae.ts.t[-1], T_TRIP + T_POST, t_pre, t_post))
    if ss.TDS.busted:
        print("  SOLVER FAILED: %s" % ss.TDS.err_msg)
        print("  the trip is at t = %.4f, so the columns below hold no post-trip state"
              % T_TRIP)

    b_pre = np.array(barq["v_pre"])
    b_post = np.array(barq["v_post"])
    a_pre = np.array([v_pre[col[b]] for b in buses])
    a_post = np.array([v_post[col[b]] for b in buses])

    print("\n%-6s %8s %8s  %8s %8s %10s"
          % ("bus", "Barq pre", "ANDES", "Barq post", "ANDES", "diff"))
    print("-" * 56)
    for k in np.argsort(-np.abs(b_post - b_pre))[:8]:
        print("%-6d %8.4f %8.4f  %8.4f %8.4f %10.4f"
              % (buses[k], b_pre[k], a_pre[k], b_post[k], a_post[k], a_post[k] - b_post[k]))
    print("-" * 56)
    print("max |pre-trip difference|  = %.4e" % np.abs(a_pre - b_pre).max())
    print("max |post-trip difference| = %.4e" % np.abs(a_post - b_post).max())
    print("min |V| post-trip: Barq %.4f, ANDES %.4f" % (b_post.min(), a_post.min()))

    t = np.array(ss.dae.ts.t)
    v_all = np.array(ss.dae.ts.y)[:, ss.Bus.v.a]
    slack_bus = int(ss.Slack.bus.v[0])
    print("trajectory over %.1f s after the trip: min |V| = %.4f, final min |V| = %.4f, "
          "slack held at %.4f"
          % (T_POST, v_all[t >= T_TRIP].min(), v_all[-1].min(), v_all[-1, col[slack_bus]]))
    return ss, v_pf, col, bus_ids


def main():
    barq = json.load(open(BARQ_RESULT))
    print("Barq reference: line %d (bus %d - bus %d), pre-trip current %.4f pu"
          % (barq["line"], barq["bus1"], barq["bus2"], barq["i_pre"]))

    compare(barq, "andes")
    ss, v_pf, col, bus_ids = compare(barq, "barq")

    # Hand the Barq-equivalent post-trip state back for a Barq-side residual check
    out = os.path.join(REPO, "results_topology/line14_andes.json")
    t = np.array(ss.dae.ts.t)
    k = int(np.argmax(t >= T_TRIP + TSTEP))
    y = np.array(ss.dae.ts.y)
    json.dump({"t": float(t[k]), "buses": bus_ids,
               "v": [float(x) for x in y[k, ss.Bus.v.a]],
               "a": [float(x) for x in y[k, ss.Bus.a.a]]},
              open(out, "w"), indent=2)
    print("\nANDES post-trip state (Barq-equivalent machines) written to %s" % out)

    # what the shunts Barq ignores are worth at the operating point
    _, v_sh = run_andes(build_case(include_shunt=True, machine_base="barq"))
    print("min power-flow |V| without shunts %.4f, with shunts %.4f, max shift %.4f"
          % (v_pf.min(), v_sh.min(), np.abs(v_sh - v_pf).max()))


if __name__ == "__main__":
    main()
