# run after init_dynamic_sim.jl

using OrdinaryDiffEq
using LinearAlgebra


n_gens = length(models.generator.bus)
n_lines = length(models.line.idx)
n_faults = length(models.fault.bus)
n_buses = length(models.bus.idx)
# TODO: add loads but dont add a load if its connected to a slack bus
n_loads = length(models.load.bus)

## state variables
idx_delta = 1:n_gens
idx_omega = idx_delta[end]+1 : idx_delta[end]+n_gens
# idx_gen_id = idx_omega[end]+1 : idx_omega[end]+n_gens
# idx_gen_iq = idx_gen_id[end]+1 : idx_gen_id[end]+n_gens
idx_line_id = idx_omega[end]+1 : idx_omega[end]+n_lines
idx_line_iq = idx_line_id[end]+1 : idx_line_id[end]+n_lines
idx_fault_id = idx_line_iq[end]+1 : idx_line_iq[end]+n_faults
idx_fault_iq = idx_fault_id[end]+1 : idx_fault_id[end]+n_faults
# idx_gen_id = idx_fault_iq[end]+1 : idx_fault_iq[end]+n_gens
# idx_gen_iq = idx_gen_id[end]+1 : idx_gen_id[end]+n_gens
idx_balance_d = idx_fault_iq[end]+1 : idx_fault_iq[end]+(length(non_slack_buses))
idx_balance_q = idx_balance_d[end]+1 : idx_balance_d[end]+(length(non_slack_buses))
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
    line_id = u[address["line_id"]]
    line_iq = u[address["line_iq"]]

    bus_vd = Vector{T}(bus.vd)  # converts Float64 -> T safely
    bus_vq = Vector{T}(bus.vq)

    bus_vd[non_slack_buses] .= u[address["balance_d"]]
    bus_vq[non_slack_buses] .= u[address["balance_q"]]

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

    address, models, _ = p

    bus = models.bus
    line = models.line

    line_id = u[address["line_id"]]
    line_iq = u[address["line_iq"]]
    

    bus_vd = Vector{T}(bus.vd)  # converts Float64 -> T safely
    bus_vq = Vector{T}(bus.vq)

    bus_vd[non_slack_buses] .= u[address["balance_d"]]
    bus_vq[non_slack_buses] .= u[address["balance_q"]]

    du[address["line_id"]] = @. bus_vd[line.bus1_idx] - bus_vd[line.bus2_idx] - line.R * line_id + line.X * line_iq
    du[address["line_iq"]] = @. bus_vq[line.bus1_idx] - bus_vq[line.bus2_idx] - line.R * line_iq - line.X * line_id

end

function solve_fault!(du, u, p)
    T = eltype(u)
    Ω = T(2*pi*60)

    address, models, _ = p

    bus = models.bus
    line = models.line
    fault = models.fault

    fault_id = u[address["fault_id"]]
    fault_iq = u[address["fault_iq"]]

    Ls = @. T(fault.l_s) / Ω

    bus_vd = Vector{T}(bus.vd)  # converts Float64 -> T safely
    bus_vq = Vector{T}(bus.vq)

    bus_vd[non_slack_buses] .= u[address["balance_d"]]
    bus_vq[non_slack_buses] .= u[address["balance_q"]]

    # du[address["fault_id"]] = @. bus_vd[fault.bus] - fault.r_s * fault_id + fault.x_s * fault_iq
    # du[address["fault_iq"]] = @. bus_vq[fault.bus] - fault.r_s * fault_iq - fault.x_s * fault_id
    du[address["fault_id"]] = @. bus_vd[fault.bus] - fault.r_s * fault_id
    du[address["fault_iq"]] = @. bus_vq[fault.bus] - fault.r_s * fault_iq

end

function balance!(du, u, p)
    T = eltype(u)
    address, models, incidence_matrix = p

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
    fault_id = u[address["fault_id"]]
    fault_iq = u[address["fault_iq"]]
    balance_d = u[address["balance_d"]]
    balance_q = u[address["balance_q"]]

    bus_vd[non_slack_buses] .= balance_d
    bus_vq[non_slack_buses] .= balance_q

    id = zeros(T, length(bus.idx))
    iq = zeros(T, length(bus.idx))

    ## load current = conjugate(load power) / conjugate(bus voltage)
    i_load = @. (load.p - 1im * load.q) / (bus_vd[load.bus] - 1im * bus_vq[load.bus])

    id[generator.bus] += @. gen_id * cos(gen_delta - pi/2) - gen_iq * sin(gen_delta - pi/2)
    id[load.bus] -= @. real(i_load)
    id[fault.bus] -= @. fault_id
    id[:] += incidence_matrix * line_id

    iq[generator.bus] += @. gen_id * sin(gen_delta - pi/2) + gen_iq * cos(gen_delta - pi/2)
    iq[load.bus] -= @. imag(i_load)
    iq[fault.bus] -= @. fault_iq
    iq[:] += incidence_matrix * line_iq

    B_mat = build_B_matrix(models)
    i_cap_d = zeros(T, length(bus.idx))
    i_cap_q = zeros(T, length(bus.idx))

    id[:] += B_mat * bus_vq
    iq[:] -= B_mat * bus_vd

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
p = (address, models, incidence_matrix)

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

    ##vd 
    i = 1
    for idx in collect(address["balance_d"])
        mass_matrix[idx, idx] = C_eq[i]
        i += 1
    end

    ##vq bus balances
    i = 1
    for idx in collect(address["balance_q"])
        mass_matrix[idx, idx] = C_eq[i]
        i += 1
    end

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

    # d#elta
    for i in collect(address["delta"])
        mass_matrix[i,i] = 1.0
    end
    # o#mega
    idx = 1
    for i in collect(address["omega"])
        mass_matrix[i,i] = models.generator.M[idx]
        idx += 1
    end

    # ## fault id
    # i = 1
    # for idx in collect(address["fault_id"])
    #     mass_matrix[idx, idx] = models.fault.l_s[i]
    #     i += 1
    # end

    # # #fault iq
    # i = 1
    # for idx in collect(address["fault_iq"])
    #     mass_matrix[idx, idx] = models.fault.l_s[i]
    #     i += 1
    # end
    
end

tspan = (0.0, 1.0)
prob0 = ODEFunction(solve_dynamic_sim!, mass_matrix=mass_matrix)
prob = ODEProblem(prob0, u0, tspan, p)

sol = solve(prob, Trapezoid(), adaptive=true, dt = 50e-5, reltol=1e-6, abstol=1e-6)

# dt = 50e-5
# prob = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, tspan, p, mass_matrix)
# sol = MyDiffEq.Solve(prob, dt, method=:Euler, adaptive=false);

using Plots

plot(sol, idxs = address["delta"])
plot(sol, idxs = address["omega"])
plot(sol, idxs = address["gen_id"])
plot(sol, idxs = address["gen_iq"])
plot(sol, idxs = address["line_id"])
plot(sol, idxs = address["line_iq"])
plot(sol, idxs = address["fault_id"])
plot(sol, idxs = address["fault_iq"])
plot(sol, idxs = address["balance_d"])
plot(sol, idxs = address["balance_q"])

fault_d = [sol_u[address["fault_id"]][1] for sol_u in sol.u]


plot(fault_d)

begin
    u1 = [sol_u[1] for sol_u in sol.u]
    u2 = [sol_u[2] for sol_u in sol.u]
    u3 = [sol_u[3] for sol_u in sol.u]
    u4 = [sol_u[4] for sol_u in sol.u]
    u5 = [sol_u[5] for sol_u in sol.u]
    u6 = [sol_u[6] for sol_u in sol.u]
    u7 = [sol_u[7] for sol_u in sol.u]
    u8 = [sol_u[8] for sol_u in sol.u]
    u9 = [sol_u[9] for sol_u in sol.u]
    u10 = [sol_u[10] for sol_u in sol.u]
    u11 = [sol_u[11] for sol_u in sol.u]
    u12 = [sol_u[12] for sol_u in sol.u]
    u13 = [sol_u[13] for sol_u in sol.u]
    u14 = [sol_u[14] for sol_u in sol.u]
    u49 = [sol_u[49] for sol_u in sol.u]
end

plot(u1)
plot!(u2)
plot!(u3)
plot!(u4)
plot(u5)
plot(u6)
plot(u7)
plot(u8)
plot(u9)
plot(u10)
plot(u11)
plot(u12)
plot(u13)
plot(u14)




omega_idx = address["omega"][1]
omega_pre_fault = [sol_u[omega_idx] for sol_u in sol.u]

models.fault.r_s = [1e5]
models.fault.x_s = [1e5]
models.fault.l_s = models.fault.x_s ./ (2*pi*60)

# fault id
i = 1
for idx in collect(address["fault_id"])
    mass_matrix[idx, idx] = models.fault.l_s[i]
    i += 1
end

# fault iq
i = 1
for idx in collect(address["fault_iq"])
    mass_matrix[idx, idx] = models.fault.l_s[i]
    i += 1
end


u0 = sol.u[end]
tspan = (0.0, 0.1)
prob0 = ODEFunction(solve_dynamic_sim!, mass_matrix=mass_matrix)
prob = ODEProblem(prob0, u0, tspan, p)

sol = solve(prob, Rodas5(), adaptive=true, dt = 50e-6, abstol=1e-4)
# sol = solve(prob, Trapezoid(), adaptive=true, dt = 50e-5, reltol=1e-6)

# dt = 50e-8
# prob = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, tspan, p, mass_matrix)
# sol = MyDiffEq.Solve(prob, dt, method=:Euler, adaptive=false);

omega_fault = [sol_u[omega_idx] for sol_u in sol.u]

plot(sol, idxs = address["delta"])
plot(sol, idxs = address["omega"])
plot(sol, idxs = address["gen_id"])
plot(sol, idxs = address["gen_iq"])
plot(sol, idxs = address["line_id"])
plot(sol, idxs = address["line_iq"])
plot(sol, idxs = address["fault_id"])
plot(sol, idxs = address["fault_iq"])
plot(sol, idxs = address["balance_d"])
plot(sol, idxs = address["balance_q"])

models.fault.r_s = [1e10]
# models.fault.x_s = [1e10]
# models.fault.l_s = models.fault.x_s ./ (2*pi*60)

# fault id
# i = 1
# for idx in collect(address["fault_id"])
#     mass_matrix[idx, idx] = models.fault.l_s[i]
#     i += 1
# end

# # fault iq
# i = 1
# for idx in collect(address["fault_iq"])
#     mass_matrix[idx, idx] = models.fault.l_s[i]
#     i += 1
# end

u0 = sol.u[end]
tspan = (0.0, 10.0)
# prob0 = ODEFunction(solve_dynamic_sim!, mass_matrix=mass_matrix)
# prob = ODEProblem(prob0, u0, tspan, p)
# sol = solve(prob, Trapezoid(), adaptive=false, dt = 50e-5)

dt = 50e-6
prob = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, tspan, p, mass_matrix)
sol = MyDiffEq.Solve(prob, dt, method=:Euler, adaptive=false);

omega_post_fault = [sol_u[omega_idx] for sol_u in sol.u]

omega = vcat(omega_pre_fault, omega_fault, omega_post_fault)
plot(omega[1:140000])

using Plots
plot(sol, idxs = address["delta"])
plot(sol, idxs = address["omega"])
plot(sol, idxs = address["gen_id"])
plot(sol, idxs = address["gen_iq"])
plot(sol, idxs = address["line_id"])
plot(sol, idxs = address["line_iq"])
plot(sol, idxs = address["fault_id"])
plot(sol, idxs = address["fault_iq"])
plot(sol, idxs = address["balance_d"])
plot(sol, idxs = address["balance_q"])