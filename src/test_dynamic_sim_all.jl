# run after init_dynamic_sim.jl
using OrdinaryDiffEq
using LinearAlgebra


n_gens = length(models.generator.bus)
n_lines = length(models.line.idx)
n_buses = length(models.bus.idx)
n_loads = length(models.load.bus)
n_faults = length(models.fault.bus)

## state variables
idx_delta = 1:n_gens
idx_omega = idx_delta[end]+1 : idx_delta[end]+n_gens
idx_line_id = idx_omega[end]+1 : idx_omega[end]+n_lines
idx_line_iq = idx_line_id[end]+1 : idx_line_id[end]+n_lines
idx_fault_id = idx_line_iq[end]+1 : idx_line_iq[end]+n_faults
idx_fault_iq = idx_fault_id[end]+1 : idx_fault_id[end]+n_faults
idx_balance_d = idx_fault_iq[end]+1 : idx_fault_iq[end]+length(non_slack_buses)
idx_balance_q = idx_balance_d[end]+1 : idx_balance_d[end]+length(non_slack_buses)
idx_gen_id = idx_balance_q[end]+1 : idx_balance_q[end]+n_gens
idx_gen_iq = idx_gen_id[end]+1 : idx_gen_id[end]+n_gens



address = Dict(
    "delta" => idx_delta,
    "omega" => idx_omega,
    "gen_id" => idx_gen_id,
    "gen_iq" => idx_gen_iq,
    "line_id" => idx_line_id,
    "line_iq" => idx_line_iq,
    "fault_id" => idx_fault_id,
    "fault_iq" => idx_fault_iq,
    "balance_d" => idx_balance_d,
    "balance_q" => idx_balance_q
)


function solve_generator!(du, u, p)
    T = eltype(u)

    Ω = T(2*pi*60)
    d = T(1.0)

    address, models, _ , _ = p

    bus = models.bus
    generator = models.generator
    line = models.line
    load = models.load
    slack = models.slack

    gen_delta = u[address["delta"]]
    gen_omega = u[address["omega"]]
    gen_id = u[address["gen_id"]]
    gen_iq = u[address["gen_iq"]]
    bus_vd = Vector{T}(bus.vd)
    bus_vq = Vector{T}(bus.vq)
    bus_vd[non_slack_buses] = @. u[address["balance_d"]]
    bus_vq[non_slack_buses] = @. u[address["balance_q"]]

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

function solve_line!(du, u, p)
    T = eltype(u)
    Ω = T(2*pi*60)

    address, models, _ , _ = p

    bus = models.bus
    line = models.line

    line_id = u[address["line_id"]]
    line_iq = u[address["line_iq"]]
    

    bus_vd = Vector{T}(bus.vd)  # converts Float64 -> T safely
    bus_vq = Vector{T}(bus.vq)

    bus_vd[non_slack_buses] = @. u[address["balance_d"]]
    bus_vq[non_slack_buses] = @. u[address["balance_q"]]

    du[address["line_id"]] = @. bus_vd[line.bus1_idx] - bus_vd[line.bus2_idx] - line.R * line_id + line.X * line_iq
    du[address["line_iq"]] = @. bus_vq[line.bus1_idx] - bus_vq[line.bus2_idx] - line.R * line_iq - line.X * line_id

end

function solve_fault!(du, u, p)
    T = eltype(u)
    Ω = T(2*pi*60)

    address, models, _, _ = p

    bus = models.bus
    line = models.line
    fault = models.fault

    fault_id = u[address["fault_id"]]
    fault_iq = u[address["fault_iq"]]

    Ls = @. T(fault.l_s) / Ω

    bus_vd = Vector{T}(bus.vd)  # converts Float64 -> T safely
    bus_vq = Vector{T}(bus.vq)

    bus_vd[non_slack_buses] = @. u[address["balance_d"]]
    bus_vq[non_slack_buses] = @. u[address["balance_q"]]

    du[address["fault_id"]] = @. bus_vd[fault.bus] - fault.r_s * fault_id + fault.x_s * fault_iq
    du[address["fault_iq"]] = @. bus_vq[fault.bus] - fault.r_s * fault_iq - fault.x_s * fault_id
    # du[address["fault_id"]] = @. bus_vd[fault.bus] - fault.r_s * fault_id
    # du[address["fault_iq"]] = @. bus_vq[fault.bus] - fault.r_s * fault_iq

end

function balance!(du, u, p)
    T = eltype(u)
    address, models, incidence_matrix, C_eq = p

    bus = models.bus
    line = models.line
    generator = models.generator
    fault = models.fault
    load = models.load
    slack = models.slack

    bus_vd = Vector{T}(bus.vd)  # converts Float64 -> T safely
    bus_vq = Vector{T}(bus.vq)

    gen_delta = u[address["delta"]]
    gen_omega = u[address["omega"]]
    gen_id = u[address["gen_id"]]
    gen_iq = u[address["gen_iq"]]
    line_id = u[address["line_id"]]
    line_iq = u[address["line_iq"]]

    bus_vd[non_slack_buses] = @. u[address["balance_d"]]
    bus_vq[non_slack_buses] = @. u[address["balance_q"]]

    fault_id = u[address["fault_id"]]
    fault_iq = u[address["fault_iq"]]

    id = zeros(T, length(bus.idx))
    iq = zeros(T, length(bus.idx))

    ## load current = conjugate(load power) / conjugate(bus voltage)
    i_load = @. (load.p - 1im * load.q) / (bus_vd[load.bus] - 1im * bus_vq[load.bus])

    id[generator.bus] += @. gen_id * sin(gen_delta) + gen_iq * cos(gen_delta)
    id[load.bus] -= @. real(i_load)
    id[fault.bus] -= @. fault_id
    id[:] += incidence_matrix * line_id

    iq[generator.bus] += @. gen_iq * sin(gen_delta) - gen_id * cos(gen_delta)
    iq[load.bus] -= @. imag(i_load)
    iq[fault.bus] -= @. fault_iq
    iq[:] += incidence_matrix * line_iq

    omega = 2*pi*60
    w1 = @. omega*C_eq * bus_vq
    w2 = @. omega*C_eq * bus_vd
    id[:] += w1
    iq[:] -= w2

    du[address["balance_d"]] = @. id[non_slack_buses]
    du[address["balance_q"]] = @. iq[non_slack_buses]
end


function solve_dynamic_sim!(du, u, p, t)
    solve_generator!(du, u, p)
    solve_line!(du, u, p)
    solve_fault!(du, u, p)
    balance!(du, u, p)
end


incidence_matrix = build_incidence_matrix(models)


du = zeros(4*n_gens + 2*n_lines + 2*n_faults + 2*length(non_slack_buses))

begin
    u0 = zeros(length(du))
    u0[address["delta"]] = models.generator.delta
    u0[address["omega"]] = models.generator.omega
    u0[address["gen_id"]] = models.generator.i_d
    u0[address["gen_iq"]] = models.generator.i_q
    u0[address["line_id"]] = models.line.i_d
    u0[address["line_iq"]] = models.line.i_q
    u0[address["fault_id"]] = [0.0]
    u0[address["fault_iq"]] = [0.0]
    u0[address["balance_d"]] = models.bus.vd[non_slack_buses]
    u0[address["balance_q"]] = models.bus.vq[non_slack_buses]
end

mass_matrix = zeros(length(du), length(du))

b_mat = build_B_matrix(models)
C_mat = b_mat ./ (2*pi*60)
C_eq = diag(C_mat)


begin

    for (i, bus_idx) in enumerate(address["balance_d"])
        # address["balance_d"] contains the indices in the state vector u
        # We need the ACTUAL bus index to get the correct C_eq
        real_bus_idx = non_slack_buses[i] 
        
        # Assign Mass Matrix
        # Vd Row
        mass_matrix[bus_idx, bus_idx] = C_eq[real_bus_idx]
        
        # Vq Row (the corresponding q index)
        q_idx = address["balance_q"][i]
        mass_matrix[q_idx, q_idx] = C_eq[real_bus_idx]
    end

    # #line id
    i = 1
    for idx in collect(address["line_id"])
        mass_matrix[idx, idx] = models.line.L[i]
        # mass_matrix[idx, idx] = 0.0
        i += 1
    end

    # #line iq
    i = 1
    for idx in collect(address["line_iq"])
        mass_matrix[idx, idx] = models.line.L[i]
        # mass_matrix[idx, idx] = 0.0
        i += 1
    end

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

p = (address, models, incidence_matrix, C_eq)
tspan = (0.0, 10.0)
prob0 = ODEFunction(solve_dynamic_sim!, mass_matrix=mass_matrix)
prob = ODEProblem(prob0, u0, tspan, p)

# sol = solve(prob, Trapezoid(), adaptive=false, dt = 50e-5, reltol=1e-6, abstol=1e-6)
# sol = solve(prob, ImplicitEuler(nlsolve = NLNewton(always_new=true)), adaptive=false, dt = 50e-5)
# sol = solve(prob, Trapezoid(), adaptive=true, dt = 50e-5, reltol=1e-6, abstol=1e-6)
sol = solve(prob, Trapezoid(), adaptive=false, dt = 50e-3)

using Plots
plot(sol, idxs = address["delta"])
plot(sol, idxs = address["omega"])
plot(sol, idxs = address["gen_id"])
plot(sol, idxs = address["gen_iq"])
plot(sol, idxs = address["line_id"])
plot(sol, idxs = address["line_iq"])
plot(sol, idxs = address["balance_d"])
plot(sol, idxs = address["balance_q"])

v1d = [sol_u[address["balance_d"]][1] for sol_u in sol.u]
v1q = [sol_u[address["balance_q"]][1] for sol_u in sol.u]

v2d = [sol_u[address["balance_d"]][2] for sol_u in sol.u]
v2q = [sol_u[address["balance_q"]][2] for sol_u in sol.u]

v3d = [sol_u[address["balance_d"]][3] for sol_u in sol.u]
v3q = [sol_u[address["balance_q"]][3] for sol_u in sol.u]

v4d = [sol_u[address["balance_d"]][4] for sol_u in sol.u]
v4q = [sol_u[address["balance_q"]][4] for sol_u in sol.u]

v5d = [sol_u[address["balance_d"]][5] for sol_u in sol.u]
v5q = [sol_u[address["balance_q"]][5] for sol_u in sol.u]

v6d = [sol_u[address["balance_d"]][6] for sol_u in sol.u]
v6q = [sol_u[address["balance_q"]][6] for sol_u in sol.u]

v7d = [sol_u[address["balance_d"]][7] for sol_u in sol.u]
v7q = [sol_u[address["balance_q"]][7] for sol_u in sol.u]

v8d = [sol_u[address["balance_d"]][8] for sol_u in sol.u]
v8q = [sol_u[address["balance_q"]][8] for sol_u in sol.u]

v9d = [sol_u[address["balance_d"]][9] for sol_u in sol.u]
v9q = [sol_u[address["balance_q"]][9] for sol_u in sol.u]

v10d = [sol_u[address["balance_d"]][10] for sol_u in sol.u]
v10q = [sol_u[address["balance_q"]][10] for sol_u in sol.u]

v11d = [sol_u[address["balance_d"]][11] for sol_u in sol.u]
v11q = [sol_u[address["balance_q"]][11] for sol_u in sol.u]

v12d = [sol_u[address["balance_d"]][12] for sol_u in sol.u]
v12q = [sol_u[address["balance_q"]][12] for sol_u in sol.u]

v13d = [sol_u[address["balance_d"]][13] for sol_u in sol.u]
v13q = [sol_u[address["balance_q"]][13] for sol_u in sol.u]

plot(abs.(v1d + 1im * v1q))
plot(abs.(v2d + 1im * v2q))
plot(abs.(v3d + 1im * v3q))
plot(abs.(v4d + 1im * v4q))
plot(abs.(v5d + 1im * v5q))
plot(abs.(v6d + 1im * v6q))
plot(abs.(v7d + 1im * v7q))
plot(abs.(v8d + 1im * v8q))
plot(abs.(v9d + 1im * v9q))
plot(abs.(v10d + 1im * v10q))
plot(abs.(v11d + 1im * v11q))
plot(abs.(v12d + 1im * v12q))
plot(abs.(v13d + 1im * v13q))
