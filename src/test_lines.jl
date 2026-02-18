# run after init_dynamic_sim.jl

using OrdinaryDiffEq
using LinearAlgebra


# n_gens = length(models.generator.bus)
n_lines = length(models.line.idx)
# n_faults = length(models.fault.bus)

## state variables
idx_line_id = 1:n_lines
idx_line_iq = n_lines+1 : 2*n_lines



address = Dict(
    "line_id" => idx_line_id,
    "line_iq" => idx_line_iq,
)


function solve_line!(du, u, p)
    T = eltype(u)
    Ω = T(2*pi*60)

    address, models, _ = p

    bus = models.bus
    line = models.line

    line_id = u[address["line_id"]]
    line_iq = u[address["line_iq"]]
    

    bus_vd = Vector{T}(bus.vd)  # converts Float64 -> T safely
    bus_vq = Vector{T}(bus.vq)

    bus_vd[non_slack_buses] = @. models.bus.vd[non_slack_buses]
    bus_vq[non_slack_buses] = @. models.bus.vq[non_slack_buses]

    du[address["line_id"]] = @. bus_vd[line.bus1_idx] - bus_vd[line.bus2_idx] - line.R * line_id + line.X * line_iq
    du[address["line_iq"]] = @. bus_vq[line.bus1_idx] - bus_vq[line.bus2_idx] - line.R * line_iq - line.X * line_id

end


function solve_dynamic_sim!(du, u, p, t)
    solve_line!(du, u, p)
end


incidence_matrix = build_incidence_matrix(models)
p = (address, models, incidence_matrix)

du = zeros(2*n_lines)

begin
    u0 = zeros(length(du))
    u0[address["line_id"]] = models.line.i_d
    u0[address["line_iq"]] = models.line.i_q
end

mass_matrix = zeros(length(du), length(du))

begin

    # #line id
    i = 1
    for idx in collect(address["line_id"])
        mass_matrix[idx, idx] = models.line.L[i]
        i += 1
    end

    # #line iq
    i = 1
    for idx in collect(address["line_iq"])
        mass_matrix[idx, idx] = models.line.L[i]
        i += 1
    end
    
end

tspan = (0.0, 1.0)
prob0 = ODEFunction(solve_dynamic_sim!, mass_matrix=mass_matrix)
prob = ODEProblem(prob0, u0, tspan, p)

# sol = solve(prob, Trapezoid(), adaptive=true, dt = 50e-5, reltol=1e-6, abstol=1e-6)
sol = solve(prob, Rodas5(), adaptive=false, dt = 50e-5)

# dt = 50e-5
# prob = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, tspan, p, mass_matrix)
# sol = MyDiffEq.Solve(prob, dt, method=:Euler, adaptive=false);

using Plots

u1 = [sol_u[1] for sol_u in sol.u]
u2 = [sol_u[2] for sol_u in sol.u]

plot(u1)
plot(u2)
plot(sol, idxs = address["line_id"])
plot(sol, idxs = address["line_iq"])



