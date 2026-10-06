# Three two-panel figures for the letter, from the CSVs written by scripts/trace_divergence.jl
# (run with h=0.001 T=0.3 tag=fig).
#   Fig. 1 (case 1, 39-bus, Newton -> LV root): (a) time-domain |V_f| pre-fault, fault-on, post-clearing
#          for direct integration / NR re-init (LV) and HC / flat re-init (HV); (b) continuation path.
#   Fig. 2 (case 2, 39-bus, Newton stalls): (a) Newton residuals: direct integration at three h, NR re-init,
#          first DAE step after HC re-init; (b) time-domain |V_f| from the HC and flat-start roots.
#   Fig. 3 (case 3, 118-bus, only HC recovers): (a) time-domain |V_f| for all methods; (b) continuation path.
#
#   julia --project=. scripts/plot_letter_figs.jl [case1=<dir>] [case2=<dir>] [case3=<dir>] [tpost=0.1] [out=figures/review1/letter]
using Plots, LaTeXStrings, Printf
args = Dict(String(split(a, "="; limit=2)[1]) => String(split(a, "="; limit=2)[2]) for a in ARGS if occursin("=", a))
D1 = get(args, "case1", "outputs/clearing_lowv/divergence_fault_bus_6_xf_0_002_clear_0_1_lfbase_power_fig")
D2 = get(args, "case2", "outputs/clearing_lowv/divergence_fault_bus_16_xf_0_001_cleared_by_trip_16_21_lfbase_power_fig")
D3 = get(args, "case3", "outputs/clearing_lowv/divergence_fault_bus_48_xf_0_001_cleared_by_trip_48_49_lfbase_power_fig")
TPOST = parse(Float64, get(args, "tpost", "0.5"))          # post-clearing window shown in the time-domain panels (s)
OUT = get(args, "out", "figures/review1/letter")
NIT = parse(Int, get(args, "nit", "10"))

function readcsv(dir, name)
    lines = filter(!isempty, readlines(joinpath(dir, name)))
    cols = Dict(String(strip(c)) => i for (i, c) in enumerate(split(lines[1], ",")))
    data = reduce(vcat, [permutedims(parse.(Float64, split(l, ","))) for l in lines[2:end]])
    (; data, cols)
end
col(t, name) = Float64.(t.data[:, t.cols[name]])
has(dir, name) = isfile(joinpath(dir, name))

default(fontfamily="Computer Modern", linewidth=3.5, markersize=8, markerstrokewidth=1.0, markerstrokecolor=:black,
        legendfontsize=14, guidefontsize=15, tickfontsize=16, titlefontsize=15, dpi=300, framestyle=:box, grid=false)
C_BLUE, C_ORANGE, C_GREEN, C_SKY, C_PURPLE, C_BLACK, C_GREY = "#0072B2", "#D55E00", "#009E73", "#56B4E9", "#CC79A7", "#333333", "#888888"
hlabel(h) = h >= 1e-3 ? @sprintf("%g ms", 1e3h) : @sprintf("%g μs", 1e6h)
ms(t) = 1e3 .* t

# ---------------------------------------------------------------- building blocks
# pre-fault + fault-on segment (t <= 0), then post-clearing runs (t >= 0), time in ms
function time_panel(dir, runs; title="(a)", legendpos=(0.75, 0.2))
    fo = readcsv(dir, "faulton_run.csv"); v_clear = col(fo, "V_faultbus")[end]     # inherited value at t = 0
    p = plot(ms(col(fo, "t")), col(fo, "V_faultbus"), color=C_BLACK, label="", xlabel="Time (ms)",
             ylabel=L"|V_f|\ \mathrm{(p.u.)}", title=title, ylims=(-0.05, 1.1), legend=legendpos)
    vline!(p, [0.0], color=C_BLACK, linestyle=:dot, linewidth=1.5, label="")
    for (file, lab, c, ls, lw) in runs
        has(dir, file) || continue
        r = readcsv(dir, file); t = col(r, "t"); v = col(r, "V_faultbus"); sel = t .<= TPOST
        # start every post-clearing curve at the inherited value, so the re-initialization
        # jump at t = 0 is drawn as a vertical segment connecting the two trajectories
        plot!(p, ms(vcat(0.0, t[sel])), vcat(v_clear, v[sel]), color=c, linestyle=ls, linewidth=lw, label=lab)
    end
    plot!(p, xlims=(-150, 1e3*TPOST))
    # annotate!(p, [(-100, 1.1, text("pre-fault", 11, :center)), (-50, 1.1, text("fault-on", 11, :center)), (0.5e3*TPOST, 1.1, text("post-clearing", 11, :center))])
    p
end
function continuation_panel(dir; title="(b)")
    hc = readcsv(dir, "hc_path.csv"); λ = col(hc, "lambda"); v = col(hc, "V_faultbus")
    p = plot(λ, v, color=C_BLUE, label="", xlabel="Continuation parameter (λ)", ylabel=L"|V_f|\ \mathrm{(p.u.)}",
             title=title, ylims=(-0.09, 1.1), xticks=0:0.25:1)
    sel = unique(round.(Int, range(1, length(λ), length=12)))
    scatter!(p, λ[sel], v[sel], color=C_BLUE, label="")
    p
end
function residual_panel(dir; title="(a)")
    p = plot(xlabel="NR iteration", ylabel=L"\|R(\mathbf{z})\|_2", yscale=:log10, title=title, legend=(0.61,0.45), legendfontsize=14)
    dfiles = sort(filter(f -> startswith(f, "direct_h") && endswith(f, ".csv"), readdir(dir)), by = f -> -parse(Float64, f[9:end-4]))
    styles = [(C_BLUE, :circle, :solid), (C_ORANGE, :diamond, :dash), (C_PURPLE, :utriangle, :dashdot)]
    kmax = 1
    for (i, file) in enumerate(dfiles)
        h = parse(Float64, file[9:end-4]); c, m, ls = styles[mod1(i, 3)]
        r = col(readcsv(dir, file), "G_2"); k = 1:min(NIT, length(r)); kmax = max(kmax, length(k))
        plot!(p, k, r[k], color=c, marker=m, linestyle=ls, label="No re-init. ($(hlabel(h)))")
    end
    r = col(readcsv(dir, "newton_reinit.csv"), "g_2"); k = 1:min(NIT, length(r)); kmax = max(kmax, length(k))
    plot!(p, k, r[k], color=C_BLACK, marker=:square, linestyle=:dashdotdot, label="NR re-init.", alpha=0.8)
    # first DAE step after the HC and the flat-start re-initialization (largest h available)
    for (prefix, lab, c, m, ls, msz) in (("post_hc_step_h", "HC re-init.", C_GREEN, :dtriangle, :dashdotdot, 7),
                                        ("post_flat_step_h", "1st step after flat re-init.", C_SKY, :circle, :dot, 5))
        post = filter(f -> startswith(f, prefix) && endswith(f, ".csv"), readdir(dir))
        isempty(post) && continue
        file = first(sort(post, by = f -> -parse(Float64, f[length(prefix)+1:end-4])))
        h = parse(Float64, file[length(prefix)+1:end-4]); r = col(readcsv(dir, file), "G_2")
        plot!(p, 1:length(r), r, color=c, marker=m, markersize=msz, linestyle=ls, label="$(lab) ($(hlabel(h)))")
        # break
    end
    hline!(p, [1e-6], color=C_BLACK, linestyle=:dot, linewidth=1.5, label="Tolerance")
    plot!(p, xlims=(0.5, kmax + 0.5))
    # plot!(xticks=([1,2,4,6,8,10]))
    # plot!()
    p
end

mkpath(dirname(OUT * "_x"))
sz = (1300, 430); mg = (left_margin=9Plots.mm, right_margin=5Plots.mm, bottom_margin=12Plots.mm, top_margin=3Plots.mm)

# ---------------------------------------------------------------- Fig. 1: silent convergence to the LV root
if isdir(D1)
    runs1 = [("direct_run.csv", "No re-init.", C_BLACK, :solid, 3.0),
             ("lv_run.csv", "NR re-init.", C_ORANGE, :dash, 3.0),
             ("hv_run.csv", "HC re-init.", C_BLUE, :solid, 3.0),
             ("flat_run.csv", "Flat-start re-init.", C_SKY, :dashdot, 2.5)]
    f1 = plot(time_panel(D1, runs1), continuation_panel(D1), layout=(1, 2), size=sz; mg...)
    # savefig(f1, OUT * "_fig1_lvroot_39.pdf"); savefig(f1, OUT * "_fig1_lvroot_39.png"); println("saved ", OUT * "_fig1_lvroot_39.pdf")
end
# ---------------------------------------------------------------- Fig. 2: Newton stall
if isdir(D2)
    runs2 = [("hv_run.csv", "HC re-init.", C_BLUE, :solid, 3.0),
             ("flat_run.csv", "Flat re-init.", C_ORANGE, :dashdot, 2.5)]
    f2 = plot(residual_panel(D2), time_panel(D2, runs2; title="(b)"), layout=(1, 2), size=sz; mg...)
    # savefig(f2, OUT * "_fig2_stall_39.pdf"); savefig(f2, OUT * "_fig2_stall_39.png"); println("saved ", OUT * "_fig2_stall_39.pdf")
end
# ---------------------------------------------------------------- Fig. 3: only HC recovers
if isdir(D3)
    runs3 = [("direct_run.csv", "No re-init.", C_PURPLE, :dashdotdot, 3.0),
             ("lv_run.csv", "NR re-init.", C_ORANGE, :dash, 3.0),
             ("flat_run.csv", "Flat-start re-init.", C_SKY, :dashdot, 2.5),
             ("hv_run.csv", "HC re-init.", C_BLUE, :solid, 3.0)]
    f3 = plot(time_panel(D3, runs3; legendpos=(0.65, 0.5)), continuation_panel(D3), layout=(1, 2), size=sz; mg...)
    # savefig(f3, OUT * "_fig3_onlyhc_118.pdf"); savefig(f3, OUT * "_fig3_onlyhc_118.png"); println("saved ", OUT * "_fig3_onlyhc_118.pdf")
end

plot(residual_panel(D3))