# this file is used to update power_flow_traditional.jl with improved methods

using Pkg
Pkg.activate(".")

include("Models.jl")
using .Models

include("data_loader.jl")
using .DataLoader

include("utils.jl")
using .Utils

# push!(LOAD_PATH, "/Users/aali27/Work/Research/DAE_Reinitialization/DAEReinitialization/MyDiffEq/MyDiffEq")
# using MyDiffEq
# data_file = "cases/SMIB_Chow/SMIB_RL_Line_DrCui.xlsx"
# data_file = "cases/SMIB_Chow/ieee14_simplified.xlsx"
data_file = "cases/Fault_cases/ieee14_fault_barq.xlsx"
# data_file = "cases/Fault_cases/ieee39_fault.xlsx"
# data_file = "cases/Fault_cases/SMIB_RL_Line_DrCui.xlsx"
# data_file = "cases/Simple_cases/ieee14_fault_barq.xlsx"

models = load_data(data_file)
n_bus = length(models.bus.idx)

# create Y matrix
Y = build_Y_matrix(models)

# get bus types
# v_update_buses = setdiff(models.bus.idx, models.generator.bus, models.slack.bus)
non_slack_buses = setdiff(models.bus.idx, models.slack.bus)
v_update_buses = setdiff(non_slack_buses, models.generator.bus)

V = copy(models.bus.v)
theta = copy(models.bus.theta)


function power_flow_traditional!(du, u, p)
    T = eltype(u)

    address, models, Y = p

    bus = models.bus
    generator = models.generator
    load = models.load

    V = Vector{T}(bus.v)
    theta = Vector{T}(bus.theta)

    V[v_update_buses] = @. u[address["balance_reactive_power"]]
    theta[non_slack_buses] = @. u[address["balance_real_power"]]

    P = zeros(T, length(bus.idx))
    Q = zeros(T, length(bus.idx))

    v = V .* exp.(1im .* theta)
    S = v .* conj.(Y * v)

    P = real.(S)
    Q = imag.(S)

    ## power balance : (power injection into network) + (load) - (generation) = 0 
    P[generator.bus] -= @. generator.p_m
    P[load.bus] += @. load.p
    Q[load.bus] += @. load.q

    du[address["balance_reactive_power"]] .= @. Q[v_update_buses]
    du[address["balance_real_power"]] .= @. P[non_slack_buses]
end



idx_balance_reactive_power = 1 : length(v_update_buses)
# activate power balance eqautions are required at both PV and PQ buses
idx_balance_real_power = idx_balance_reactive_power[end]+1 : idx_balance_reactive_power[end]+(length(non_slack_buses))

# 2 vars per line + 2 vars at PQ buses + 1 var at PV buses
# du = zeros(2*(length(v_update_buses)) + length(models.generator.bus))
address = Dict(
    "balance_reactive_power" => idx_balance_reactive_power,
    "balance_real_power" => idx_balance_real_power
)

# u0 = [1.0044, 0.31948, 0.18429]
# u0 = [1.0, 0.0, 0.0]
u0 = vcat(models.bus.v[v_update_buses], models.bus.theta[non_slack_buses])
du = zeros(length(u0))

using NonlinearSolve

function powerflow!(u, p)
    T = eltype(u)
    du = zeros(T, length(u))
    power_flow_traditional!(du, u, p)
    return du
end


p = (address, models, Y)

prob = NonlinearProblem(powerflow!, u0, p)
sol = solve(prob, abstol = 1e-9, reltol = 1e-9)

# populate models with the solution
models.bus.v[v_update_buses] .= sol.u[address["balance_reactive_power"]]
models.bus.theta[non_slack_buses] .= sol.u[address["balance_real_power"]]

@show models.bus.v
@show models.bus.theta


function phasor2DP!(bus)
    vdq = @. bus.v * exp(1im * bus.theta)
    bus.vd = real(vdq)
    bus.vq = imag(vdq)
end
# convert to d-q components
phasor2DP!(models.bus)

# compute reactive power at buses
v_complex = @. models.bus.v * exp(1im * models.bus.theta)
s_bus = v_complex .* conj(Y * v_complex)
p_bus = real(s_bus)
q_bus = imag(s_bus)

# power balance at buses
#  power_gen - power_load - power_line = 0 => power_gen = power_load + power_line
p_gen = zeros(length(models.bus.idx))
p_gen[models.load.bus] += @. models.load.p
p_gen += @. p_bus
@show p_gen[models.generator.bus]

q_gen = zeros(length(models.bus.idx))
q_gen[models.load.bus] += @. models.load.q
q_gen += @. q_bus

@show q_gen[models.generator.bus]