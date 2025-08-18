using Pkg
Pkg.activate(".")

include("Models.jl")
using .Models

include("data_loader.jl")
using .DataLoader

include("powerflow.jl")
using .PowerFlow

system, sol = run_powerflow_example()

system.models.fault.l_s = 1e10 .* (system.models.fault.l_s ./ (2*pi*60))
system.models.fault.r_s =  system.models.fault.r_s .+ 1e10 


n_gen = length(system.models.generator.bus)
n_line = length(system.models.line.idx)
n_bus = length(system.models.bus.idx)
n_fault = length(system.models.fault.bus)
n_load = length(system.models.load.bus)

du = zeros(5*n_gen)

idx_delta = 1:n_gen
idx_omega = idx_delta[end]+1 : idx_delta[end]+n_gen
idx_stator_d = idx_omega[end]+1 : idx_omega[end]+n_gen
idx_stator_q = idx_stator_d[end]+1 : idx_stator_d[end]+n_gen
idx_e_q_prime = idx_stator_q[end]+1 : idx_stator_q[end]+ n_gen

address = Dict(
    "delta" => idx_delta,
    "omega" => idx_omega,
    "stator_d" => idx_stator_d,
    "stator_q" => idx_stator_q,
    "e_q_prime" => idx_e_q_prime,
)


# system.models.bus.vd[1] = 0.9778807983015356
# system.models.bus.vq[1] = 0.3234951998301541
# system.models.bus.vd[3] = 0.9874149369646668
# system.models.bus.vq[3] = 0.1840576136964672
# system.models.line.i_d[1] = 0.9140869456592352
# system.models.line.i_q[1] = 0.15496961898679878
# system.models.line.i_d[2] = 0.9140869456592352
# system.models.line.i_q[2] = 0.15496961898679878

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


    gen_delta = u[address["delta"]]
    gen_omega = u[address["omega"]]
    gen_id = u[address["stator_d"]]
    gen_iq = u[address["stator_q"]]
    gen_eq_prime = u[address["e_q_prime"]]


    # make working copies with the right eltype
    bus_vd = Vector{T}(bus.vd)  # converts Float64 -> T safely
    bus_vq = Vector{T}(bus.vq)



    du[address["delta"]] = @. Ω * (gen_omega - one(T))
    du[address["omega"]] = @. (generator.p_m - (gen_id * bus_vd[generator.bus] * sin(gen_delta) - 
                                     gen_id * bus_vq[generator.bus] * cos(gen_delta) + 
                                     gen_iq * bus_vd[generator.bus] * cos(gen_delta) + 
                                     gen_iq * bus_vq[generator.bus] * sin(gen_delta)) - 
                                     d * (gen_omega - one(T))) / (generator.M)
end


function stator_equation!(du, u, address::Dict, models :: NamedTuple)
    T = eltype(u)

    bus = models.bus
    line = models.line
    generator = models.generator
    fault = models.fault
    load = models.load
    slack = models.slack



    gen_delta = u[address["delta"]]
    gen_omega = u[address["omega"]]
    gen_id = u[address["stator_d"]]
    gen_iq = u[address["stator_q"]]
    gen_eq_prime = u[address["e_q_prime"]]

    # make working copies with the right eltype
    bus_vd = Vector{T}(bus.vd)  # converts Float64 -> T safely
    bus_vq = Vector{T}(bus.vq)




    du[address["stator_d"]] = @. gen_eq_prime - gen_id * generator.x_d_prime - bus_vd[generator.bus] * cos(gen_delta) - bus_vq[generator.bus] * sin(gen_delta)
    du[address["stator_q"]] = @. gen_iq * generator.x_d_prime - bus_vd[generator.bus] * sin(gen_delta) + bus_vq[generator.bus] * cos(gen_delta)
end



function balance_equation!(du, u, address::Dict, models :: NamedTuple)
    T = eltype(u)

    bus = models.bus
    line = models.line
    generator = models.generator
    fault = models.fault
    load = models.load
    slack = models.slack


    gen_delta = u[address["delta"]]
    gen_omega = u[address["omega"]]
    gen_id = u[address["stator_d"]]
    gen_iq = u[address["stator_q"]]
    gen_eq_prime = u[address["e_q_prime"]]

    # make working copies with the right eltype
    bus_vd = Vector{T}(bus.vd)  # converts Float64 -> T safely
    bus_vq = Vector{T}(bus.vq)






    real_power = zeros(T, length(bus.idx))



    p_load = @. load.p

    # p_fault = @. bus_vd[fault.bus] * fault_id + bus_vq[fault.bus] * fault_iq

    p_gen = @. gen_id * bus_vd[generator.bus] * sin(gen_delta) - 
    gen_id * bus_vq[generator.bus] * cos(gen_delta) + 
    gen_iq * bus_vd[generator.bus] * cos(gen_delta) + 
    gen_iq * bus_vq[generator.bus] * sin(gen_delta)

    # p_gen = @. generator.p_m


    p_line_from = @. bus_vd[line.bus1_idx] * line.i_d + bus_vq[line.bus1_idx] * line.i_q
    p_line_to = @. bus_vd[line.bus2_idx] * line.i_d + bus_vq[line.bus2_idx] * line.i_q
   

    real_power[load.bus] -= p_load
    # real_power[fault.bus] -= p_fault
    real_power[generator.bus] += p_gen
    real_power[line.bus1_idx] -= p_line_from
    real_power[line.bus2_idx] += p_line_to


    du[address["e_q_prime"]] = @. real_power[generator.bus]
end


using NonlinearSolve

u0 = ones(5*n_gen)


function initialize_dynamic_models!(u, p)
    T = eltype(u)
    du = zeros(T, 5*n_gen)
    address, models = p

    # states = unpack_u(u)
    
    swing_equation!(du, u, address, models)
    stator_equation!(du, u, address, models)
    balance_equation!(du, u, address, models)

    return du
end

p = (address, system.models)

prob = NonlinearProblem(initialize_dynamic_models!, u0, p)
sol = solve(prob, reltol = 1e-10)


