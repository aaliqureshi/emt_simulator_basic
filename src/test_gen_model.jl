# test the generator dynamic model
using OrdinaryDiffEq
using LinearAlgebra


n_gens = length(models.generator.bus)
# n_lines = length(models.line.idx)
# n_faults = length(models.fault.bus)
# n_buses = length(models.bus.idx)
# TODO: add loads but dont add a load if its connected to a slack bus
# n_loads = length(models.load.bus)

## state variables
idx_delta = 1:n_gens
idx_omega = idx_delta[end]+1 : idx_delta[end]+n_gens
idx_gen_id = idx_omega[end]+1 : idx_omega[end]+n_gens
idx_gen_iq = idx_gen_id[end]+1 : idx_gen_id[end]+n_gens

address = Dict(
    "delta" => idx_delta,
    "omega" => idx_omega,
    "gen_id" => idx_gen_id,
    "gen_iq" => idx_gen_iq,
)

function solve_generator!(du, u, p)
    T = eltype(u)

    Ω = T(2*pi*60)
    d = T(1.0)

    address, models, _ = p

    bus = models.bus
    generator = models.generator
    line = models.line
    load = models.load
    slack = models.slack

    gen_delta = u[address["delta"]]
    gen_omega = u[address["omega"]]
    gen_id = u[address["gen_id"]]
    gen_iq = u[address["gen_iq"]]

    bus_vd = Vector{T}(bus.vd)  # converts Float64 -> T safely
    bus_vq = Vector{T}(bus.vq)

    bus_vd[non_slack_buses] = @. models.bus.vd[non_slack_buses]
    bus_vq[non_slack_buses] = @. models.bus.vq[non_slack_buses]

    du[address["delta"]] = @. Ω * (gen_omega - one(T))
    du[address["omega"]] = @. generator.p_m - (gen_id * bus_vd[generator.bus] * sin(gen_delta) - 
                                     gen_id * bus_vq[generator.bus] * cos(gen_delta) + 
                                     gen_iq * bus_vd[generator.bus] * cos(gen_delta) + 
                                     gen_iq * bus_vq[generator.bus] * sin(gen_delta)) - 
                                     d * (gen_omega - one(T))
    du[address["gen_id"]] = @. generator.e_q_prime - gen_id * generator.x_d_prime - 
                                bus_vd[generator.bus] * cos(gen_delta) - 
                                bus_vq[generator.bus] * sin(gen_delta)
    du[address["gen_iq"]] = @. gen_iq * generator.x_d_prime - 
                                bus_vd[generator.bus] * sin(gen_delta) + 
                                bus_vq[generator.bus] * cos(gen_delta)
    
end

function solve_dynamic_sim!(du, u, p, t)
    solve_generator!(du, u, p)
    # solve_line!(du, u, p)
    # solve_fault!(du, u, p)
    # balance!(du, u, p)
end

du = zeros(4*n_gens)
incidence_matrix = build_incidence_matrix(models)
p = (address, models, incidence_matrix)

begin
    u0 = zeros(length(du))
    u0[address["delta"]] = models.generator.delta
    u0[address["omega"]] = models.generator.omega
    u0[address["gen_id"]] = models.generator.i_d
    u0[address["gen_iq"]] = models.generator.i_q
end

mass_matrix = zeros(length(u0), length(u0))

begin
    # delta
    for i in collect(address["delta"])
        mass_matrix[i,i] = 1.0
    end
    # omega
    idx = 1
    for i in collect(address["omega"])
        mass_matrix[i,i] = models.generator.M[idx]
        idx += 1
    end
    
end

tspan = (0.0, 1.0)
prob0 = ODEFunction(solve_dynamic_sim!, mass_matrix=mass_matrix)
prob = ODEProblem(prob0, u0, tspan, p)

# sol = solve(prob, Trapezoid(), adaptive=true, dt = 50e-5, reltol=1e-6, abstol=1e-6)
sol = solve(prob, Rodas5(), adaptive=false, dt = 50e-5)

u1 = [sol_u[1] for sol_u in sol.u]

using Plots

plot(sol, idxs = address["delta"])
plot(sol, idxs = address["omega"])
plot(sol, idxs = address["gen_id"])
plot(sol, idxs = address["gen_iq"])

plot(u1)