using Plots, CairoMakie, CSV, DataFrames

function read_correction_norm_data()
    cases = ["SMIB", "14bus", "39bus", "wecc"]
    dts = [0.0005, 5e-5, 5e-6]
    data = Dict{String, Dict{Float64, Vector{Float64}}}()
    for case in cases
        data[case] = Dict{Float64, Vector{Float64}}()
        for dt in dts
            if case == "wecc" && dt == 5e-6
                continue
            end
            file_name = "$(case)_$(dt)_correction_norm.csv"
            data[case][dt] = CSV.read(file_name, DataFrame)[:, :correction_norm]
        end
    end
    return data
end

data = read_correction_norm_data()


function plotter(data, g)
    dts = [0.0005, 5e-5, 5e-6]
    n1 = data[0.0005][2:g]
    n2 = data[5e-5][2:g]
    if haskey(data, 5e-6)
        n3 = data[5e-6][2:g]
    end
    k = collect(2:g)

    fig = Figure(
        size = (950, 650),
        fontsize = 22,
        backgroundcolor = :white
    )

    ax = Axis(
        fig[1, 1],
        xlabel = "Newton iteration",
        ylabel = "Newton correction norm",
        yscale = log10,
        xlabelsize = 24,
        ylabelsize = 24,
        xticklabelsize = 18,
        yticklabelsize = 18,
        xgridvisible = true,
        ygridvisible = true,
        spinewidth = 1.5
    )

    if haskey(data, 5e-6)
        lines!(ax, k, n3, linewidth = 3, label = "Δt = 5×10^-6")
    end
    # scatter!(ax, k1, n1, markersize = 9)

    lines!(ax, k, n2, linewidth = 3, label = "Δt = 5×10^-5")
    # scatter!(ax, k2, n2, markersize = 9)

    lines!(ax, k, n1, linewidth = 3, label = "Δt = 5×10^-4")
    # scatter!(ax, k3, n3, markersize = 9)

    axislegend(
        ax,
        position = :lt,
        framevisible = true,
        padding = (8, 8, 8, 8),
        labelsize = 18
    )

    resize_to_layout!(fig)

    return fig
end


# for case in ["SMIB", "14bus", "39bus", "wecc"]
#     fig_name = "newton_correction_norm_$(case).pdf"
#     g = 10
#     fig = plotter(data[case], g)
#     save(fig_name, fig)
#     save(fig_name, fig, px_per_unit = 3)
# end

fig = plotter(data["14bus"], 10)
save("newton_correction_norm_14bus.pdf", fig)
