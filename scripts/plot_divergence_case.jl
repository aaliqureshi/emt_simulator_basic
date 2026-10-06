# Figure for the letter: Newton non-convergence at fault clearing on the 39-bus system
# (ZIP loads without PQBRAK, power-mismatch bus equations, RMS network), from the
# CSVs written by scripts/trace_divergence.jl.
#   (a) continuation path |V_f|(lambda) of the homotopy re-initialization
#   (b) Newton residual vs iteration: direct integration (first step, three h),
#       Newton re-init from the inherited state, flat-start re-init
#   (c) |V_f| along the Newton iterates (re-init) and along the homotopy path
#
#   julia --project=. scripts/plot_divergence_case.jl [dir=outputs/clearing_lowv/divergence_fault_bus_16_xf_0_001_cleared_by_trip_16_21_lfbase] [out=figures/review1/divergence_39_clearing.pdf] [nit=30]
using Plots, LaTeXStrings, Printf
args = Dict(String(split(a, "="; limit=2)[1]) => String(split(a, "="; limit=2)[2]) for a in ARGS if occursin("=", a))
dir = get(args, "dir", "outputs/clearing_lowv/divergence_fault_bus_16_xf_0_001_cleared_by_trip_16_21_lfbase_power")
out = get(args, "out", "figures/review1/divergence_39_clearing.pdf")
NIT = parse(Int, get(args, "nit", "30"))
fault_bus = match(r"bus_(\d+)", dir).captures[1]

function readcsv(name)
    lines = filter(!isempty, readlines(joinpath(dir, name)))
    cols = Dict(String(strip(c)) => i for (i, c) in enumerate(split(lines[1], ",")))
    data = reduce(vcat, [permutedims(parse.(Float64, split(l, ","))) for l in lines[2:end]])
    (; data, cols)
end
col(t, name) = Float64.(t.data[:, t.cols[name]])

default(fontfamily="Computer Modern", linewidth=3.0, markersize=7, markerstrokewidth=1.0, markerstrokecolor=:black,
        legendfontsize=12, guidefontsize=15, tickfontsize=14, titlefontsize=15, dpi=300, framestyle=:box, grid=false)
C_BLUE, C_ORANGE, C_GREEN, C_SKY, C_PURPLE, C_BLACK = "#0072B2", "#D55E00", "#009E73", "#56B4E9", "#CC79A7", "#333333"

# (a) continuation path
hc = readcsv("hc_path.csv")
λ, vhc = col(hc, "lambda"), col(hc, "V_faultbus")
pa = plot(λ, vhc, color=C_BLUE, label="", xlabel="Continuation parameter (λ)", ylabel=L"|V_f|\ \mathrm{(p.u.)}",
          title="(a) Continuation path", ylims=(0, 1.15), xticks=0:0.25:1)
sel = unique(round.(Int, range(1, length(λ), length=12)))
scatter!(pa, λ[sel], vhc[sel], color=C_BLUE, label="")

# (b) residual histories
pb = plot(xlabel="NR iteration", ylabel=L"\|R(\mathbf{z})\|_2", yscale=:log10, title="(b) Newton residual", legend=:bottomright)
kmax = 1
# direct-integration first-step histories: one per direct_h*.csv, labelled by h
dfiles = sort(filter(f -> startswith(f, "direct_h") && endswith(f, ".csv"), readdir(dir)),
              by = f -> -parse(Float64, f[9:end-4]))
styles = [(C_BLUE, :circle, :solid), (C_ORANGE, :diamond, :dash), (C_PURPLE, :utriangle, :dashdot)]
hlabel(h) = h >= 1e-3 ? @sprintf("%g ms", 1e3h) : @sprintf("%g μs", 1e6h)
for (i, file) in enumerate(dfiles)
    h = parse(Float64, file[9:end-4]); c, m, ls = styles[mod1(i, length(styles))]
    t = readcsv(file); r = col(t, "G_2"); k = 1:min(NIT, length(r))
    plot!(pb, k, r[k], color=c, marker=m, linestyle=ls, label="No re-init. (h = $(hlabel(h)))"); global kmax = max(kmax, length(k))
end
nr = readcsv("newton_reinit.csv"); r = col(nr, "g_2"); k = 1:min(NIT, length(r))
plot!(pb, k, r[k], color=C_BLACK, marker=:square, linestyle=:solid, label="NR re-init."); kmax = max(kmax, length(k))
if isfile(joinpath(dir, "flat_reinit.csv"))
    fl = readcsv("flat_reinit.csv"); r = col(fl, "g_inf")
    plot!(pb, 1:length(r), r, color=C_SKY, marker=:dtriangle, linestyle=:dashdotdot, label="Flat-start re-init."); kmax = max(kmax, length(r))
end
hline!(pb, [1e-6], color=C_BLACK, linestyle=:dot, linewidth=1.5, label="Tolerance")
plot!(pb, xlims=(0.5, kmax + 0.5))

# (c) if Newton converged to a different root (lv_run.csv exists): time-domain |V_f| from
#     Newton's root and from the homotopy root; otherwise |V_f| along the Newton iterates
if isfile(joinpath(dir, "lv_run.csv"))
    lv = readcsv("lv_run.csv"); hv = readcsv("hv_run.csv")
    pc = plot(col(hv, "t"), col(hv, "V_faultbus"), color=C_BLUE, label="From homotopy root", xlabel="Time (s)",
              ylabel=L"|V_f|\ \mathrm{(p.u.)}", title="(c) Post-clearing trajectory", legend=:right, ylims=(0, 1.15))
    plot!(pc, col(lv, "t"), col(lv, "V_faultbus"), color=C_BLACK, linestyle=:dash, label="From NR / flat-start root")
else
    v_nr = col(nr, "V_faultbus"); k = 1:min(NIT, length(v_nr))
    pc = plot(k, v_nr[k], color=C_BLACK, marker=:square, label="NR re-init. iterates", xlabel="NR iteration",
              ylabel=L"|V_f|\ \mathrm{(p.u.)}", title="(c) Faulted-bus voltage", legend=:topright, ylims=(0, max(1.15, 1.1*maximum(v_nr[k]))), xlims=(0.5, NIT + 0.5))
    hline!(pc, [vhc[end]], color=C_BLUE, linestyle=:dash, linewidth=2, label="Homotopy root")
    hline!(pc, [vhc[1]], color=C_ORANGE, linestyle=:dot, linewidth=2, label="Inherited (fault-on) value")
end

plt = plot(pa, pb, pc, layout=(1, 3), size=(1900, 480), left_margin=10Plots.mm, right_margin=3Plots.mm, bottom_margin=13Plots.mm, top_margin=3Plots.mm)
mkpath(dirname(out)); savefig(plt, out)
plt2 = plot(pa, pb, layout=(1, 2), size=(1300, 420), left_margin=10Plots.mm, right_margin=3Plots.mm, bottom_margin=13Plots.mm, top_margin=3Plots.mm)
savefig(plt2, replace(out, ".pdf" => "_2panel.pdf"))
savefig(plt, replace(out, ".pdf" => ".png"))
println("saved ", out, " and the 2-panel / png variants")
