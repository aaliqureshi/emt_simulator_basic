# TODO: add provision for different ZIP coefficients for P and Q

module BusModel

export Bus, balance!, phasor2DP!, motor_d_load

mutable struct Bus{T<:Real}
    idx :: Vector{Int32}
    orig_idx :: Vector{Int32}
    v :: Vector{T}
    theta :: Vector{T}
    vd:: Vector{T}
    vq:: Vector{T}

    function Bus{T}(n::Integer) where {T<:Real}
        new{T}(Vector{Int32}(undef, n),
               Vector{Int32}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n))
    end
end

function phasor2DP!(bus)
    vdq = @. bus.v * exp(1im * bus.theta)
    bus.vd = real(vdq)
    bus.vq = imag(vdq)
end

#=
# Previous balance! implementation, preserved for comparison.
function balance!(du, u, p)
    T = eltype(u)
    address, models, incidence_matrix, C_eq, non_slack_buses, _ = p

    bus = models.bus
    generator = models.generator
    fault = models.fault
    load = models.load

    bus_vd = Vector{T}(bus.vd)
    bus_vq = Vector{T}(bus.vq)

    gen_delta = u[address["delta"]]
    gen_id = u[address["gen_id"]]
    gen_iq = u[address["gen_iq"]]
    line_id = u[address["line_id"]]
    line_iq = u[address["line_iq"]]

    bus_vd[non_slack_buses] = @. u[address["balance_d"]]
    bus_vq[non_slack_buses] = @. u[address["balance_q"]]

    fault_id = u[address["fault_id"]]
    fault_iq = u[address["fault_iq"]]

    id = zeros(T, length(bus.idx))
    iq = zeros(T, length(bus.idx))

    # k_pz = zeros(T, length(load.bus))
    # k_pi = zeros(T, length(load.bus))
    # k_pp = zeros(T, length(load.bus))

    # k_pz[:] .= 0.7
    # k_pi[:] .= 0.1
    # k_pp[:] .= 0.2
    k_pz = T(0.7)
    k_pi = T(0.1)
    k_pp = T(0.2)

    vd_l = bus_vd[load.bus]
    vq_l = bus_vq[load.bus]

    # used only by the current-injection model below
    # v_min = T(0.7)
    # p_l  = Vector{T}(load.p)
    # q_l  = Vector{T}(load.q)

    # initial operating point
    v0_pq = models.bus.v[models.load.bus]

    # current point
    v_pq = @. hypot(vd_l, vq_l)

    # kp_fac = ones(T, length(load.bus))
    # v_min = 0.7

    # for idx in eachindex(v_pq)
    #     if v_pq[idx] < v_min
    #         kp_fac[idx] = 0.4881 - 0.4999*cos(3.964*v_pq[idx]) + 0.1389*sin(3.964*v_pq[idx])
    #     elseif v_pq[idx] >= v_min
    #         kp_fac[idx] = 1.0 
    #     end
    # end

    # k_pp = @. k_pp * kp_fac

    # Z component: I_z = conj(S₀) · V / |V₀|²
    # i_d_z = @. k_pz * (p_l * vd_l + q_l * vq_l) / v0_sq
    # i_q_z = @. k_pz * (p_l * vq_l - q_l * vd_l) / v0_sq

    # I component: I_i = conj(S₀) · V / (|V₀| · |V|)
    # v_sq_i = @. max(v_sq, v_min^2)
    # denom_i = @. sqrt(v0_sq * v_sq_i)
    # i_d_i = @. k_pi * (p_l * vd_l + q_l * vq_l) / denom_i
    # i_q_i = @. k_pi * (p_l * vq_l - q_l * vd_l) / denom_i

    # P component: I_p = conj(S₀) · V / |V|²
    # v_sq_p = @. max(v_sq, v_min^2)
    # i_d_p = @. k_pp * (p_l * vd_l + q_l * vq_l) / v_sq_p
    # i_q_p = @. k_pp * (p_l * vq_l - q_l * vd_l) / v_sq_p

    ## ZIP Milano model pg. 315
    p0 = models.load.p
    q0 = models.load.q

    p_load = @. p0 * (
        k_pz * (v_pq / v0_pq)^2 + 
        k_pi * (v_pq / v0_pq) +
        k_pp
        )

    q_load = @. q0 * (
        k_pz * (v_pq / v0_pq)^2 +
        k_pi * (v_pq / v0_pq) +
        k_pp
        )

    # Total load current (only needed by the current-injection formulation above)
    # i_load_d = @. i_d_z + i_d_i + i_d_p
    # i_load_q = @. i_q_z + i_q_i + i_q_p

    id[generator.bus] += @. gen_id * sin(gen_delta) + gen_iq * cos(gen_delta)
    # id[load.bus] -= i_load_d
    id[fault.bus] -= @. fault_id
    id[:] += incidence_matrix * line_id

    iq[generator.bus] += @. gen_iq * sin(gen_delta) - gen_id * cos(gen_delta)
    # iq[load.bus] -= i_load_q
    iq[fault.bus] -= @. fault_iq
    iq[:] += incidence_matrix * line_iq

    omega = 2*pi*60
    omega_C = @. omega * C_eq
    w1 = @. omega_C * bus_vq
    w2 = @. omega_C * bus_vd
    id[:] += w1
    iq[:] -= w2

    p_h = zeros(T, length(bus.idx))
    q_h = zeros(T, length(bus.idx))

    p_h[non_slack_buses] = @. id[non_slack_buses] * bus_vd[non_slack_buses] + iq[non_slack_buses] * bus_vq[non_slack_buses]
    q_h[non_slack_buses] = @. id[non_slack_buses] * bus_vq[non_slack_buses] - iq[non_slack_buses] * bus_vd[non_slack_buses]
    p_h[models.load.bus] .-= p_load
    q_h[models.load.bus] .-= q_load

    # du[address["balance_d"]] = @. id[non_slack_buses]
    # du[address["balance_q"]] = @. iq[non_slack_buses]

    du[address["balance_d"]] = @. p_h[non_slack_buses]
    du[address["balance_q"]] = @. q_h[non_slack_buses]
end

# =#

#
"""
    _openipsl_load_factors(v, pqbrak, characteristic)

Return the constant-power and constant-current multipliers from OpenIPSL's
Electrical.Loads.PSSE.BaseClasses.baseLoad. Voltage and PQBRAK are absolute p.u.
Characteristic 1 is piecewise quadratic; characteristic 2 uses the trigonometric
power reduction and exponential current reduction from the upstream model.

The strict inequalities of characteristic 1 are reproduced intentionally:
at exactly v = 0 or v = PQBRAK/2, upstream falls through to kP = 1.
Source: https://github.com/OpenIPSL/OpenIPSL/blob/master/OpenIPSL/Electrical/Loads/PSSE/BaseClasses/baseLoad.mo
"""
function _openipsl_load_factors(v, pqbrak, characteristic)
    if characteristic == 1
        kP = if 0 < v < pqbrak / 2
            2 * (v / pqbrak)^2
        elseif pqbrak / 2 < v < pqbrak
            1 - 2 * ((v - pqbrak) / pqbrak)^2
        else
            one(v)
        end
        return kP, one(v)
    elseif characteristic == 2
        kP = v < pqbrak ? 0.4881 - 0.4999*cos(3.964*v) + 0.1389*sin(3.964*v) : one(v)
        kI = v < 0.5 ? 1.502*1.769*v^(1.769 - 1)*exp(-1.502*v^1.769) : one(v)
        return kP, kI
    end
    throw(ArgumentError("OpenIPSL characteristic must be 1 or 2"))
end

"""
WECC composite load model motor D (single-phase A/C compressor) constants.
Per unit on the motor base `S = P_init / lf`; `pf` is the compressor power factor.
Source: WECC Composite Load Model Specification (MVS, April 2021), motor D section.
"""
const MOTOR_D = (lf=1.0, pf=0.98, vbrk=0.86, vstall=0.6, rstall=0.1, xstall=0.1,
                 kp1=0.0, np1=1.0, kq1=6.0, nq1=2.0, kp2=12.0, np2=3.2, kq2=11.0, nq2=2.5)

"""
    _motor_d_vstallbrk(prm)

Voltage where the motor D run-state-II power curve meets the stall impedance
curve, searched by bisection in `[0.4, prm.vstall]`; `prm.vstall` if the curves
do not intersect there.
"""
function _motor_d_vstallbrk(prm)
    gstall = prm.rstall / (prm.rstall^2 + prm.xstall^2)
    f(v) = prm.lf + prm.kp2*(prm.vbrk - v)^prm.np2 - gstall*v^2
    lo, hi = 0.4, prm.vstall
    f(lo) * f(hi) > 0 && return prm.vstall
    for _ in 1:60
        mid = 0.5*(lo + hi)
        f(lo) * f(mid) <= 0 ? (hi = mid) : (lo = mid)
    end
    return 0.5*(lo + hi)
end

"""
    _motor_d_pq(v, v0, p_init, vstallbrk, prm)

Motor D active and reactive power (system p.u.) at voltage `v` for a motor whose
initial demand is `p_init` at voltage `v0`. Three algebraic states: run I above
`vbrk`, run II (rising demand) between `vstallbrk` and `vbrk`, stall impedance below.
"""
function _motor_d_pq(v, v0, p_init, vstallbrk, prm)
    v0 > prm.vbrk || throw(ArgumentError("motor D initialization requires v0 > vbrk"))
    S = p_init / prm.lf
    p0 = prm.lf
    q0 = p0*tan(acos(prm.pf)) - prm.kq1*(v0 - prm.vbrk)^prm.nq1
    if v > prm.vbrk
        pm = p0 + prm.kp1*(v - prm.vbrk)^prm.np1
        qm = q0 + prm.kq1*(v - prm.vbrk)^prm.nq1
    elseif v > vstallbrk
        pm = p0 + prm.kp2*(prm.vbrk - v)^prm.np2
        qm = q0 + prm.kq2*(prm.vbrk - v)^prm.nq2
    else
        z2 = prm.rstall^2 + prm.xstall^2
        pm = prm.rstall / z2 * v^2
        qm = prm.xstall / z2 * v^2
    end
    return S*pm, S*qm
end

"""
    motor_d_load(fraction; params=MOTOR_D)

Motor D load specification for the `motor` keyword of `balance!`: `fraction` of
each load's initial P is served by a WECC motor D with parameters `params`.
"""
function motor_d_load(fraction; params=MOTOR_D)
    0 <= fraction <= 1 || throw(ArgumentError("motor D fraction must be in [0, 1]"))
    return (; fraction, params, vstallbrk=_motor_d_vstallbrk(params))
end

"""
    _motor_d_load(p0, q0, v, v0, motor)

Split a load with initial demand `p0 + j q0` into its motor D share and the
remainder. Returns the motor P and Q at voltage `v`, and the remaining initial
P and Q. The motor's initial Q follows from `motor.params.pf` and is taken out of
the load's Q so the power flow point is preserved.
"""
function _motor_d_load(p0, q0, v, v0, motor)
    p_md0 = motor.fraction * p0
    q_md0 = p_md0 * tan(acos(motor.params.pf))
    p_md, q_md = _motor_d_pq(v, v0, p_md0, motor.vstallbrk, motor.params)
    return p_md, q_md, p0 - p_md0, q0 - q_md0
end

"""
    balance!(du, u, p; zip=(0.7, 0.1, 0.2), pqbrak=0.7e-6,
             characteristic=2, low_voltage=true, motor=nothing)

Power-balance bus equations with OpenIPSL PSSE low-voltage load characteristics.
`zip` gives (constant impedance, constant current, constant power) fractions of
the specified P and Q at the initial bus voltage. Set `low_voltage=false` for an
unmodified ZIP control.

OpenIPSL Load.mo uses S_Y*v^2 + kI*S_I*v + kP*S_P for both P and Q. Here those
coefficients are mapped from the existing load data as S_Y = kz*S0/v0^2,
S_I = ki*S0/v0 and S_P = kp*S0. This preserves the chosen ZIP mix rather than
applying OpenIPSL's separate default load-transfer fractions a and b.
Source: https://github.com/OpenIPSL/OpenIPSL/blob/master/OpenIPSL/Electrical/Loads/PSSE/Load.mo

`motor` adds a motor load at every load bus; `nothing` (default) means pure ZIP.
Pass `motor_d_load(fraction)` for a WECC motor D share, the remainder being ZIP.
"""
function balance!(du, u, p; zip=(0.7, 0.1, 0.2), pqbrak=0.7e-6,
                  characteristic=2, low_voltage=true, motor=nothing)
    pqbrak > 0 || throw(ArgumentError("pqbrak must be positive"))
    characteristic in (1, 2) || throw(ArgumentError("characteristic must be 1 or 2"))
    k_z, k_i, k_p = zip
    T = eltype(u)
    address, models, incidence_matrix, C_eq, non_slack_buses, _ = p
    bus = models.bus
    generator = models.generator
    fault = models.fault
    load = models.load

    bus_vd = Vector{T}(bus.vd)
    bus_vq = Vector{T}(bus.vq)
    bus_vd[non_slack_buses] = u[address["balance_d"]]
    bus_vq[non_slack_buses] = u[address["balance_q"]]

    gen_delta = u[address["delta"]]
    gen_id = u[address["gen_id"]]
    gen_iq = u[address["gen_iq"]]
    line_id = u[address["line_id"]]
    line_iq = u[address["line_iq"]]
    fault_id = u[address["fault_id"]]
    fault_iq = u[address["fault_iq"]]

    p_load = zeros(T, length(load.bus))
    q_load = zeros(T, length(load.bus))
    for idx in eachindex(load.bus)
        b = load.bus[idx]
        v = hypot(bus_vd[b], bus_vq[b])
        v0 = bus.v[b]
        p0, q0 = load.p[idx], load.q[idx]
        if motor !== nothing
            p_load[idx], q_load[idx], p0, q0 = _motor_d_load(p0, q0, v, v0, motor)
        end
        kP, kI = low_voltage ? _openipsl_load_factors(v, pqbrak, characteristic) : (one(v), one(v))
        load_factor = k_z*(v/v0)^2 + kI*k_i*(v/v0) + kP*k_p
        p_load[idx] += p0 * load_factor
        q_load[idx] += q0 * load_factor
    end

    # Preserve the existing network current assembly and power-balance residual.
    id = zeros(T, length(bus.idx))
    iq = zeros(T, length(bus.idx))
    id[generator.bus] += @. gen_id*sin(gen_delta) + gen_iq*cos(gen_delta)
    iq[generator.bus] += @. gen_iq*sin(gen_delta) - gen_id*cos(gen_delta)
    id[fault.bus] -= fault_id
    iq[fault.bus] -= fault_iq
    id[:] += incidence_matrix * line_id
    iq[:] += incidence_matrix * line_iq

    # omega = 2*pi*60
    # id[:] += @. omega*C_eq*bus_vq
    # iq[:] -= @. omega*C_eq*bus_vd

    p_h = zeros(T, length(bus.idx))
    q_h = zeros(T, length(bus.idx))
    p_h[non_slack_buses] = @. id[non_slack_buses]*bus_vd[non_slack_buses] + iq[non_slack_buses]*bus_vq[non_slack_buses]
    q_h[non_slack_buses] = @. id[non_slack_buses]*bus_vq[non_slack_buses] - iq[non_slack_buses]*bus_vd[non_slack_buses]
    p_h[load.bus] .-= p_load
    q_h[load.bus] .-= q_load
    du[address["balance_d"]] = p_h[non_slack_buses]
    du[address["balance_q"]] = q_h[non_slack_buses]

    # du[address["balance_d"]] = id[non_slack_buses]
    # du[address["balance_q"]] = iq[non_slack_buses]

    # C dv/dt = i_net - conj(S_load / V)
    # for (k, b) in enumerate(load.bus)
    #     v2 = bus_vd[b]^2 + bus_vq[b]^2
    #     id[b] -= (p_load[k]*bus_vd[b] + q_load[k]*bus_vq[b]) / v2
    #     iq[b] -= (p_load[k]*bus_vq[b] - q_load[k]*bus_vd[b]) / v2
    # end
    # du[address["balance_d"]] = id[non_slack_buses]
    # du[address["balance_q"]] = iq[non_slack_buses]
end

# =#
end # module
