using Pkg
Pkg.activate(".")

include("Models.jl")
using .Models

include("data_loader.jl")
using .DataLoader


data_file = "cases/SMIB_Chow/SMIB_RL_Line_DrCui.xlsx"

models_dict = load_data(data_file)

models = (bus=models_dict["bus"], line=models_dict["line"], generator=models_dict["generator"], fault=models_dict["fault"], load=models_dict["load"], slack=models_dict["slack"])

function phasor2DP!(bus)
    vdq = @. bus.v * exp(1im * bus.theta)
    bus.vd = real(vdq)
    bus.vq = imag(vdq)
end

phasor2DP!(models.bus)

models.line.X = models.line.X ./ (2*pi*60)


n_gen = length(models.generator.bus)
n_line = length(models.line.idx)
n_bus = length(models.bus.idx)
n_load = length(models.load.bus)
n_slack = length(models.slack.bus)

du = zeros(2*n_line + 2*(n_load) + n_gen)

vd_update_buses = sort(setdiff(models.load.bus, models.generator.bus))
vq_update_buses = sort(union(models.load.bus, models.generator.bus))

idx_line_d = 1:n_line
idx_line_q = idx_line_d[end]+1 : idx_line_d[end]+n_line
idx_balance_q = idx_line_q[end]+1 : idx_line_q[end]+(n_load)
idx_balance_d = idx_balance_q[end]+1 : idx_balance_q[end]+(n_load + n_gen)


address = Dict(
    "line_d" => idx_line_d,
    "line_q" => idx_line_q,
    "balance_d" => idx_balance_d,
    "balance_q" => idx_balance_q
)


function line_equation!(du, u, address::Dict, models :: NamedTuple)
    T = eltype(u)
    Ω = T(2*pi*60)

    bus = models.bus
    line = models.line

    line_id = u[address["line_d"]]
    line_iq = u[address["line_q"]]

    # make working copies with the right eltype
    bus_v = Vector{T}(bus.v)  # converts Float64 -> T safely
    bus_theta = Vector{T}(bus.theta)

    bus_vd = Vector{T}(bus.vd)
    bus_vq = Vector{T}(bus.vq)


    # bus_vd[vd_update_buses] .= u[address["balance_d"]]
    # bus_vq[vq_update_buses] .= u[address["balance_q"]]

    bus_v[vd_update_buses] .= u[address["balance_q"]]
    bus_theta[vq_update_buses] .= u[address["balance_d"]]

    # bus_theta[generator.bus] = @. asin(bus_vq[generator.bus] ./ bus_v[generator.bus])
    # bus_vd[generator.bus] = @. bus_v[generator.bus] * cos(bus_theta[generator.bus])

    bus_vd = @. bus_v * cos(bus_theta)
    bus_vq = @. bus_v * sin(bus_theta)


    du[address["line_d"]] = @. ((bus_vd[line.bus1_idx] - bus_vd[line.bus2_idx] - line.R * line_id) / line.X) + (Ω * line_iq)
    du[address["line_q"]] = @. ((bus_vq[line.bus1_idx] - bus_vq[line.bus2_idx] - line.R * line_iq) / line.X) - (Ω * line_id)
end

function power_balance!(du, u, address::Dict, models :: NamedTuple)
    T = eltype(u)

    bus = models.bus
    line = models.line
    generator = models.generator
    load = models.load

    line_id = u[address["line_d"]]
    line_iq = u[address["line_q"]]


    # make working copies with the right eltype
    bus_v = Vector{T}(bus.v)  # converts Float64 -> T safely
    bus_theta = Vector{T}(bus.theta)

    bus_vd = Vector{T}(bus.vd)
    bus_vq = Vector{T}(bus.vq)


    # bus_vd[vd_update_buses] .= u[address["balance_d"]]
    # bus_vq[vq_update_buses] .= u[address["balance_q"]]

    bus_v[vd_update_buses] .= u[address["balance_q"]]
    bus_theta[vq_update_buses] .= u[address["balance_d"]]

    bus_vd = @. bus_v * cos(bus_theta)
    bus_vq = @. bus_v * sin(bus_theta)

    # bus_theta[generator.bus] = @. asin(bus_vq[generator.bus] ./ bus_v[generator.bus])
    # bus_vd[generator.bus] = @. bus_v[generator.bus] * cos(bus_theta[generator.bus])


    real_power = zeros(T, length(bus.idx))
    reactive_power = zeros(T, length(bus.idx))


    p_load = load.p
    q_load = load.q


    p_gen = generator.p_m
    q_gen = generator.q_m

    p_line_from = @. bus_vd[line.bus1_idx] * line_id + bus_vq[line.bus1_idx] * line_iq
    q_line_from= @. bus_vq[line.bus1_idx] * line_id - bus_vd[line.bus1_idx] * line_iq
    p_line_to = @. bus_vd[line.bus2_idx] * line_id + bus_vq[line.bus2_idx] * line_iq
    q_line_to = @. bus_vq[line.bus2_idx] * line_id - bus_vd[line.bus2_idx] * line_iq


    real_power[load.bus] -= p_load
    # real_power[fault.bus] -= p_fault
    real_power[generator.bus] += p_gen
    real_power[line.bus1_idx] -= p_line_from
    real_power[line.bus2_idx] += p_line_to

    reactive_power[load.bus] -= q_load
    # reactive_power[fault.bus] -= q_fault
    reactive_power[generator.bus]+=q_gen
    reactive_power[line.bus1_idx] -= q_line_from
    reactive_power[line.bus2_idx] += q_line_to

    # du[address["balance_d"]] = @. real_power[vd_update_buses]
    # du[address["balance_q"]] = @. reactive_power[vq_update_buses]

    du[address["balance_d"]] = @. real_power[non_slack_buses]
    du[address["balance_q"]] = @. reactive_power[load.bus]
end



non_slack_buses = setdiff(models.bus.idx, models.slack.bus)
u0 = vcat(models.line.i_d, models.line.i_q, models.bus.v[vd_update_buses], models.bus.theta[vq_update_buses])
# u0 = vcat(models.line.i_d, models.line.i_q, models.bus.vd[vd_update_buses], models.bus.vq[vq_update_buses])
# u0 = [0.9140869456592352, 0.9140869456592352,
# 0.1549696189867898, 0.15496961898679878,
# 1.004422950206025, 0.31948009885442474, 0.1842884927989168]

du = zeros(length(u0))

using OrdinaryDiffEq


function powerflow!(du, u, p, t)
    address, models = p
    line_equation!(du, u, address, models)
    power_balance!(du, u, address, models)
end

M = zeros(length(u0), length(u0))

tspan = (0.0, 10.0)

p = (address, models)

prob0 = ODEFunction(powerflow!, mass_matrix=M)
prob = ODEProblem(prob0, u0, tspan, p)
sol = solve(prob, Trapezoid(), adaptive=false, dt = 50e-5)



@show du