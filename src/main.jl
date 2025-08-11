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
models.fault.l_s = 1e10 .* (models.fault.l_s ./ (2*pi*60))
models.fault.r_s =  models.fault.r_s .+ 1e10 

u = vcat(generator.delta, generator.omega, line.i_d, line.i_q, 
     fault.i_d, fault.i_q, generator.i_d, generator.i_q,
     bus.vd[[1,3]], bus.vq[[1,3]])


n_gen = length(models.generator.bus)
n_line = length(models.line.idx)
n_bus = length(models.bus.idx)
n_fault = length(models.fault.bus)
n_load = length(models.load.bus)

du = zeros(4*n_gen + 2*n_line + 2*n_fault + 2*(n_bus-1))

idx_swing_delta = 1:n_gen
idx_swing_omega = idx_swing_delta[end]+1 : idx_swing_delta[end]+n_gen
idx_line_d = idx_swing_omega[end]+1 : idx_swing_omega[end]+n_line
idx_line_q = idx_line_d[end]+1 : idx_line_d[end]+n_line
idx_fault_d = idx_line_q[end]+1 : idx_line_q[end]+n_fault
idx_fault_q = idx_fault_d[end]+1 : idx_fault_d[end]+n_fault
idx_stator_d = idx_fault_q[end]+1 : idx_fault_q[end]+n_gen
idx_stator_q = idx_stator_d[end]+1 : idx_stator_d[end]+n_gen
idx_balance_d = idx_stator_q[end]+1 : idx_stator_q[end]+n_bus - length(models.slack.bus)
idx_balance_q = idx_balance_d[end]+1 : idx_balance_d[end]+n_bus - length(models.slack.bus)

address = Dict(
    "swing_delta" => idx_swing_delta,
    "swing_omega" => idx_swing_omega,
    "line_d" => idx_line_d,
    "line_q" => idx_line_q,
    "fault_d" => idx_fault_d,
    "fault_q" => idx_fault_q,
    "stator_d" => idx_stator_d,
    "stator_q" => idx_stator_q,
    "balance_d" => idx_balance_d,
    "balance_q" => idx_balance_q
)

u = vcat(generator.delta, generator.omega, line.i_d, line.i_q, 
     fault.i_d, fault.i_q, generator.i_d, generator.i_q,
     bus.vd[[1,3]], bus.vq[[1,3]])



function swing_equation!(du, u, address::Dict, models :: NamedTuple)
    T = eltype(u)
    Ω = T(2*pi*60)
    d = T(1.0)

    bus = models.bus
    line = models.line
    generator = models.generator
    fault = models.fault
    load = models.load
    slack = models.slack


    gen_delta = u[address["swing_delta"]]
    gen_omega = u[address["swing_omega"]]
    gen_id = u[address["stator_d"]]
    gen_iq = u[address["stator_q"]]
    line_id = u[address["line_d"]]
    line_iq = u[address["line_q"]]
    fault_id = u[address["fault_d"]]
    fault_iq = u[address["fault_q"]]

    # make working copies with the right eltype
    bus_vd = Vector{T}(bus.vd)  # converts Float64 -> T safely
    bus_vq = Vector{T}(bus.vq)
    bd = u[address["balance_d"]]
    bq = u[address["balance_q"]]


    # Test
    non_slack_buses = setdiff(bus.idx, slack.bus)
    bus_vd[non_slack_buses] .= bd
    bus_vq[non_slack_buses] .= bq


    du[address["swing_delta"]] = @. Ω * (gen_omega - one(T))
    du[address["swing_omega"]] = @. (generator.p_m - (gen_id * bus_vd[generator.bus] * sin(gen_delta) - 
                                     gen_id * bus_vq[generator.bus] * cos(gen_delta) + 
                                     gen_iq * bus_vd[generator.bus] * cos(gen_delta) + 
                                     gen_iq * bus_vq[generator.bus] * sin(gen_delta)) - 
                                     d * (gen_omega - one(T))) / (generator.M)
end

# swing_equation!(du, u, address, models)

function line_equation!(du, u, address::Dict, models :: NamedTuple)
    T = eltype(u)
    Ω = T(2*pi*60)

    bus = models.bus
    line = models.line
    generator = models.generator
    fault = models.fault
    load = models.load
    slack = models.slack

    gen_delta = u[address["swing_delta"]]
    gen_omega = u[address["swing_omega"]]
    gen_id = u[address["stator_d"]]
    gen_iq = u[address["stator_q"]]
    line_id = u[address["line_d"]]
    line_iq = u[address["line_q"]]
    fault_id = u[address["fault_d"]]
    fault_iq = u[address["fault_q"]]

    # make working copies with the right eltype
    bus_vd = Vector{T}(bus.vd)  # converts Float64 -> T safely
    bus_vq = Vector{T}(bus.vq)
    bd = u[address["balance_d"]]
    bq = u[address["balance_q"]]


    non_slack_buses = setdiff(bus.idx, slack.bus)
    bus_vd[non_slack_buses] .= bd
    bus_vq[non_slack_buses] .= bq


    du[address["line_d"]] = @. ((bus_vd[line.bus1_idx] - bus_vd[line.bus2_idx] - line.R * line_id) / line.X) + (Ω * line_iq)
    du[address["line_q"]] = @. ((bus_vq[line.bus1_idx] - bus_vq[line.bus2_idx] - line.R * line_iq) / line.X) - (Ω * line_id)
end

# line_equation!(du, u, address, models)

function fault_equation!(du, u, address::Dict, models :: NamedTuple)
    T = eltype(u)
    Ω = T(2*pi*60)

    bus = models.bus
    line = models.line
    generator = models.generator
    fault = models.fault
    load = models.load
    slack = models.slack



    gen_delta = u[address["swing_delta"]]
    gen_omega = u[address["swing_omega"]]
    gen_id = u[address["stator_d"]]
    gen_iq = u[address["stator_q"]]
    line_id = u[address["line_d"]]
    line_iq = u[address["line_q"]]
    fault_id = u[address["fault_d"]]
    fault_iq = u[address["fault_q"]]

    # make working copies with the right eltype
    bus_vd = Vector{T}(bus.vd)  # converts Float64 -> T safely
    bus_vq = Vector{T}(bus.vq)
    bd = u[address["balance_d"]]
    bq = u[address["balance_q"]]


    non_slack_buses = setdiff(bus.idx, slack.bus)
    bus_vd[non_slack_buses] .= bd
    bus_vq[non_slack_buses] .= bq


    du[address["fault_d"]] = @. (bus_vd[fault.bus] - (fault.r_s * fault_id)) / fault.l_s + (Ω * fault_iq)
    du[address["fault_q"]] = @. (bus_vq[fault.bus] - (fault.r_s * fault_iq)) / fault.l_s - (Ω * fault_id)
end

# fault_equation!(du, u, address, models)

function stator_equation!(du, u, address::Dict, models :: NamedTuple)
    T = eltype(u)

    bus = models.bus
    line = models.line
    generator = models.generator
    fault = models.fault
    load = models.load
    slack = models.slack



    gen_delta = u[address["swing_delta"]]
    gen_omega = u[address["swing_omega"]]
    gen_id = u[address["stator_d"]]
    gen_iq = u[address["stator_q"]]
    line_id = u[address["line_d"]]
    line_iq = u[address["line_q"]]
    fault_id = u[address["fault_d"]]
    fault_iq = u[address["fault_q"]]

    # make working copies with the right eltype
    bus_vd = Vector{T}(bus.vd)  # converts Float64 -> T safely
    bus_vq = Vector{T}(bus.vq)
    bd = u[address["balance_d"]]
    bq = u[address["balance_q"]]


    non_slack_buses = setdiff(bus.idx, slack.bus)
    bus_vd[non_slack_buses] .= bd
    bus_vq[non_slack_buses] .= bq


    du[address["stator_d"]] = @. generator.e_q_prime - gen_id * generator.x_d_prime - bus_vd[generator.bus] * cos(gen_delta) - bus_vq[generator.bus] * sin(gen_delta)
    du[address["stator_q"]] = @. gen_iq * generator.x_d_prime - bus_vd[generator.bus] * sin(gen_delta) + bus_vq[generator.bus] * cos(gen_delta)
end

# stator_equation!(du, u, address, models)

function balance_equation!(du, u, address::Dict, models :: NamedTuple)
    T = eltype(u)

    bus = models.bus
    line = models.line
    generator = models.generator
    fault = models.fault
    load = models.load
    slack = models.slack


    gen_delta = u[address["swing_delta"]]
    gen_omega = u[address["swing_omega"]]
    gen_id = u[address["stator_d"]]
    gen_iq = u[address["stator_q"]]
    line_id = u[address["line_d"]]
    line_iq = u[address["line_q"]]
    fault_id = u[address["fault_d"]]
    fault_iq = u[address["fault_q"]]

    # make working copies with the right eltype
    bus_vd = Vector{T}(bus.vd)  # converts Float64 -> T safely
    bus_vq = Vector{T}(bus.vq)
    bd = u[address["balance_d"]]
    bq = u[address["balance_q"]]


    non_slack_buses = setdiff(bus.idx, slack.bus)
    bus_vd[non_slack_buses] .= bd
    bus_vq[non_slack_buses] .= bq


    real_power = zeros(T, length(bus.idx))
    reactive_power = zeros(T, length(bus.idx))

    load_id = [0.0]
    load_iq = [0.0]





    # p_load = @. bus_vd[load.bus] * load_id + bus_vq[load.bus] * load_iq
    # q_load = @. bus_vq[load.bus] * load_id - bus_vd[load.bus] * load_iq

    p_load = @. load.p
    q_load = @. load.q

    p_fault = @. bus_vd[fault.bus] * fault_id + bus_vq[fault.bus] * fault_iq
    q_fault = @. bus_vq[fault.bus] * fault_id - bus_vd[fault.bus] * fault_iq

    p_gen = @. gen_id * bus_vd[generator.bus] * sin(gen_delta) - 
    gen_id * bus_vq[generator.bus] * cos(gen_delta) + 
    gen_iq * bus_vd[generator.bus] * cos(gen_delta) + 
    gen_iq * bus_vq[generator.bus] * sin(gen_delta)

    q_gen = @. gen_id * bus_vd[generator.bus] * cos(gen_delta) + 
    gen_id * bus_vq[generator.bus] * sin(gen_delta) - 
    gen_iq * bus_vd[generator.bus] * sin(gen_delta) + 
    gen_iq * bus_vq[generator.bus] * cos(gen_delta)

    p_line_from = @. bus_vd[line.bus1_idx] * line_id + bus_vq[line.bus1_idx] * line_iq
    q_line_from= @. bus_vq[line.bus1_idx] * line_id - bus_vd[line.bus1_idx] * line_iq
    p_line_to = @. bus_vd[line.bus2_idx] * line_id + bus_vq[line.bus2_idx] * line_iq
    q_line_to = @. bus_vq[line.bus2_idx] * line_id - bus_vd[line.bus2_idx] * line_iq


    real_power[load.bus] -= p_load
    real_power[fault.bus] -= p_fault
    real_power[generator.bus] += p_gen
    real_power[line.bus1_idx] -= p_line_from
    real_power[line.bus2_idx] += p_line_to

    reactive_power[load.bus] -= q_load
    reactive_power[fault.bus] -= q_fault
    reactive_power[generator.bus]+=q_gen
    reactive_power[line.bus1_idx] -= q_line_from
    reactive_power[line.bus2_idx] += q_line_to

    du[address["balance_d"]] = @. real_power[non_slack_buses]
    du[address["balance_q"]] = @. reactive_power[non_slack_buses]
end

# balance_equation!(du, u, address, models)


function power_flow!(du, u, address::Dict, models :: NamedTuple)
    T = eltype(u)
    Ω = T(2*pi*60)
    d = T(1.0)

    bus = models.bus
    line = models.line
    generator = models.generator
    fault = models.fault
    load = models.load
    slack = models.slack


    gen_delta = u[address["swing_delta"]]
    gen_omega = u[address["swing_omega"]]
    gen_id = u[address["stator_d"]]
    gen_iq = u[address["stator_q"]]
    line_id = u[address["line_d"]]
    line_iq = u[address["line_q"]]
    fault_id = u[address["fault_d"]]
    fault_iq = u[address["fault_q"]]

    # make working copies with the right eltype
    bus_vd = Vector{T}(bus.vd)  # converts Float64 -> T safely
    bus_vq = Vector{T}(bus.vq)
    bd = u[address["balance_d"]]
    bq = u[address["balance_q"]]


    non_slack_buses = setdiff(bus.idx, slack.bus)
    bus_vd[non_slack_buses] .= bd
    bus_vq[non_slack_buses] .= bq


    real_power = zeros(T, length(bus.idx))
    reactive_power = zeros(T, length(bus.idx))

    load_id = [0.0]
    load_iq = [0.0]


    p_load = load.p
    q_load = load.q

    p_fault = zeros(length(fault.idx))
    q_fault = zeros(length(fault.idx))

    p_gen = generator.p_m
    q_gen = zeros(length(generator.idx))

    p_line_from = @. bus_vd[line.bus1_idx] * line_id + bus_vq[line.bus1_idx] * line_iq
    q_line_from= @. bus_vq[line.bus1_idx] * line_id - bus_vd[line.bus1_idx] * line_iq
    p_line_to = @. bus_vd[line.bus2_idx] * line_id + bus_vq[line.bus2_idx] * line_iq
    q_line_to = @. bus_vq[line.bus2_idx] * line_id - bus_vd[line.bus2_idx] * line_iq


    real_power[load.bus] -= p_load
    real_power[fault.bus] -= p_fault
    real_power[generator.bus] += p_gen
    real_power[line.bus1_idx] -= p_line_from
    real_power[line.bus2_idx] += p_line_to

    reactive_power[load.bus] -= q_load
    reactive_power[fault.bus] -= q_fault
    reactive_power[generator.bus]+=q_gen
    reactive_power[line.bus1_idx] -= q_line_from
    reactive_power[line.bus2_idx] += q_line_to

    du[address["balance_d"]] = @. real_power[non_slack_buses]
    du[address["balance_q"]] = @. reactive_power[non_slack_buses]
end

using OrdinaryDiffEq

# u0 = vcat(generator.delta, generator.omega, line.i_d, line.i_q, 
#      fault.i_d, fault.i_q, generator.i_d, generator.i_q,
#      bus.vd[[1,3]], bus.vq[[1,3]])

u0 = [0.5234150168964286, 1.0, 
      0.9140869456592352, 0.9140869456592352,
      0.1549696189867898, 0.15496961898679878,
      0.0, 0.0,
      0.32267613448856314, 0.8691666349239708,
      0.9778807983015356, 0.9874149369646668,
      0.3234951998301541, 0.1840576136964672]



# u0 = ones(18)      
tspan = (0.0, 10.0)

M = zeros(length(u0), length(u0))

# for i in range(1, 2*n_gen + 2*n_line)
#     M[i,i] = 1.0
# end

# M[1,1] = 1.0
# M[2,2] = 1.0

# function unpack_u(u)
#     generator.delta = u[address["swing_delta"]]
#     generator.omega = u[address["swing_omega"]]
#     generator.i_d = u[address["stator_d"]]
#     generator.i_q = u[address["stator_q"]]
# end
function fx!(du, u, p, t)
    address, models = p

    # states = unpack_u(u)
    
    swing_equation!(du, u, address, models)
    line_equation!(du, u, address, models)
    fault_equation!(du, u, address, models)
    stator_equation!(du, u, address, models)
    balance_equation!(du, u, address, models)
end

p = (address, models)

# p = ()

prob0 = ODEFunction(fx!, mass_matrix=M)
prob = ODEProblem(prob0, u0, tspan, p)
sol = solve(prob, Trapezoid(), adaptive=false, dt = 50e-5)

using Plots

plot(sol, idxs = 1)
plot(sol, idxs = 2)
plot(sol, idxs = 3)
plot(sol, idxs = 4)
plot(sol, idxs = 5)
plot(sol, idxs = 6)
plot(sol, idxs = 7)
plot(sol, idxs = 8)
plot(sol, idxs = 9)
plot(sol, idxs = 10)