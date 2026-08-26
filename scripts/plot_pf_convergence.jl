# Plot power-flow solver convergence (residual norm vs. correction norm).
#
# Usage (after running `sol = solve_power_flow!(sys)` in main.jl):
#   include("scripts/plot_pf_convergence.jl")
#   plot_pf_convergence(sol)

module SolverPlot
using Plots, LaTeXStrings

export plot_pf_convergence, plot_of_eigen_spectrum

"""
Plot the solver convergence recorded in `sol`: the residual norm on the left
y-axis and the correction norm on the right y-axis (both log10 scale),
keeping only the last `last_n` iterations. Returns the plot object.
"""



function plot_pf_convergence(sol; last_n = 50000, overlay=false, label="")
    res = sol.residual_norm
    # cor = sol.correction_norm

    # The two logs can differ in length by one, so slice each independently.
    res_idx = max(1, length(res) - last_n + 1) : length(res)

    plot_func = overlay ? plot! : plot

    # Optional: set a clean default style once at the top of your script
    default(
        fontfamily = "Computer Modern",   # or "serif" / "sans-serif"
        linewidth  = 2,
        markersize = 4,
        grid       = true,
        gridalpha  = 0.3,
        legendfontsize = 10,
        guidefontsize  = 12,
        tickfontsize   = 10,
        titlefontsize  = 12,
        size       = (500, 350),          # good aspect ratio for single-column figures
        dpi        = 300
    )

    plt = plot_func(res_idx, res[res_idx],
               label = label,
            #    color = :blue,
               xlabel = "Iteration number",
               ylabel = L"|r_k\|",
               yscale = :log10,
               legend = :topright,
               title = "Solver convergence",
               framestyle = :box,)

    return plt
end


"""
Strip plot of Hessian eigenvalues: GD (J'J) vs PGD (J'HJ).
`which` is `:init` (flat start) or `:final` (returned iterate).
"""
function plot_of_eigen_spectrum(sol; which=:init, overlay=false, label="")
    snap = which === :final ? sol.stats.final : sol.stats.init
    sol.stats.stored || error("sol.stats is empty; pass store_stats=true")

    lam_gd  = filter(>(0), snap.eigvals_JtJ)
    lam_pgd = filter(>(0), snap.eigvals_JtHJ)

    y_gd  = fill(2.0, length(lam_gd))
    y_pgd = fill(1.0, length(lam_pgd))

    plot_func = overlay ? scatter! : scatter

    plt = plot_func(lam_gd, y_gd;
                    xlabel = L"\lambda",
                    xscale = :log10,
                    ylabel = "",
                    yticks = ([1.0, 2.0], ["PGD", "GD"]),
                    ylims  = (0.4, 2.6),
                    markersize = 5,
                    markerstrokewidth = 0,
                    markeralpha = 0.7,
                    label = isempty(label) ? "GD" : label * " GD",
                    legend = :bottomright,
                    framestyle = :box,
                    title = which === :final ? "Spectrum at solution" :
                                               "Spectrum at initial guess")
    scatter!(lam_pgd, y_pgd;
             markersize = 5,
             markerstrokewidth = 0,
             markeralpha = 0.7,
             label = isempty(label) ? "PGD" : label * " PGD")
    return plt
end
end
