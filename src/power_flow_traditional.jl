using Pkg
Pkg.activate(".")

include("Models.jl")
using .Models

include("data_loader.jl")
using .DataLoader


# data_file = "cases/SMIB_Chow/SMIB_RL_Line_DrCui.xlsx"
data_file = "cases/ieee14/ieee14_simplified.xlsx"
# data_file = "cases/5bus/pjm5bus_simplified.xlsx"
# data_file = "cases/ieee39/ieee39_simplified.xlsx"

models_dict = load_data(data_file)

models = (bus=models_dict["bus"], line=models_dict["line"], generator=models_dict["generator"], fault=models_dict["fault"], load=models_dict["load"], slack=models_dict["slack"])

n_bus = length(models.bus.idx)

Y = zeros(Complex, n_bus, n_bus)

# create Y matrix conjugate
for line_idx in eachindex(models.line.idx)
    Y[models.line.bus1_idx[line_idx], models.line.bus2_idx[line_idx]] = -1 / (models.line.R[line_idx] + 1im * models.line.X[line_idx])
    Y[models.line.bus2_idx[line_idx], models.line.bus1_idx[line_idx]] = -1 / (models.line.R[line_idx] + 1im * models.line.X[line_idx])
end
for row_idx in 1:n_bus
    Y[row_idx, row_idx] = Y[row_idx, row_idx] + sum(-1 .* Y[row_idx, :])
end

Y_abs = abs.(Y)
Y_angle = angle.(Y)

println(Y_abs)
println(Y_angle)

v_update_buses = setdiff(models.bus.idx, models.generator.bus, models.slack.bus)
non_slack_buses = setdiff(models.bus.idx, models.slack.bus)

V = copy(models.bus.v)
theta = copy(models.bus.theta)


function power_flow_traditional!(du, u, address, models)
    T = eltype(u)

    bus = models.bus
    line = models.line
    generator = models.generator
    load = models.load
    slack = models.slack

    V = Vector{T}(bus.v)
    theta = Vector{T}(bus.theta)
    # theta = copy(models.bus.theta)

    V[v_update_buses] .= u[address["balance_reactive_power"]]
    theta[non_slack_buses] .= u[address["balance_real_power"]]

    P = zeros(T, length(bus.idx))
    Q = zeros(T, length(bus.idx))

    idx = 1
    for i in eachindex(bus.idx)
        for j in 1:length(V)
            P[idx] += Y_abs[i,j] * V[j] * cos(theta[i] - theta[j] - Y_angle[i,j])
            Q[idx] += Y_abs[i,j] * V[j] * sin(theta[i] - theta[j] - Y_angle[i,j])
        end
        P[idx] *= V[i]
        Q[idx] *= V[i]
        idx += 1
    end

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
du = zeros(2*(length(v_update_buses)) + length(models.generator.bus))
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
    address, models = p
    # line_equation!(du, u, address, models)
    power_flow_traditional!(du, u, address, models)
    return du
end


p = (address, models)

prob = NonlinearProblem(powerflow!, u0, p)
sol = solve(prob, abstol = 1e-9, reltol = 1e-9)

# populate models with the solution
models.bus.v[v_update_buses] .= sol.u[address["balance_reactive_power"]]
models.bus.theta[non_slack_buses] .= sol.u[address["balance_real_power"]]

@show models.bus.v
@show models.bus.theta