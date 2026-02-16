module SystemModel

export System, build_system

include("utils.jl"); using .Utils

using LinearAlgebra

mutable struct System
    models::NamedTuple
    Y::Matrix{ComplexF64}
    non_slack_buses::Vector{Int32}
    v_update_buses::Vector{Int32}
    incidence_matrix::Matrix{Int32}
    C_eq::Vector{Float64}
end

function build_system(models)
    Y = build_Y_matrix(models)
    non_slack_buses = setdiff(models.bus.idx, models.slack.bus)
    v_update_buses = setdiff(non_slack_buses, models.generator.bus)
    incidence_matrix = build_incidence_matrix(models)
    C_eq = diag(build_B_matrix(models) ./ (2*pi*60))

    return System(models, Y, non_slack_buses, v_update_buses, incidence_matrix, C_eq)
end

end # module
