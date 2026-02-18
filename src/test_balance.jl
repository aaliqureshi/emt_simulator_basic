# run after init_dynamic_sim.jl

using OrdinaryDiffEq
using LinearAlgebra


n_buses = length(models.bus.idx)

## state variables
idx_balance_d = 1 : length(non_slack_buses)
idx_balance_q = length(non_slack_buses)+1 : 2*length(non_slack_buses)



address = Dict(
    "balance_d" => idx_balance_d,
    "balance_q" => idx_balance_q
)


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

    gen_delta = models.generator.delta
    gen_omega = models.generator.omega
    gen_id = models.generator.i_d
    gen_iq = models.generator.i_q
    line_id = models.line.i_d
    line_iq = models.line.i_q

    balance_d = u[address["balance_d"]]
    balance_q = u[address["balance_q"]]

    bus_vd[non_slack_buses] .= balance_d
    bus_vq[non_slack_buses] .= balance_q

    id = zeros(T, length(bus.idx))
    iq = zeros(T, length(bus.idx))

    ## load current = conjugate(load power) / conjugate(bus voltage)
    i_load = @. (load.p - 1im * load.q) / (bus_vd[load.bus] - 1im * bus_vq[load.bus])

    # id[generator.bus] += @. gen_id * cos(gen_delta - pi/2) - gen_iq * sin(gen_delta - pi/2)
    id[generator.bus] += @. gen_id * sin(gen_delta) + gen_iq * cos(gen_delta)
    id[load.bus] -= @. real(i_load)
    # id[fault.bus] -= @. fault_id
    id[:] += incidence_matrix * line_id

    iq[generator.bus] += @. gen_iq * sin(gen_delta) - gen_id * cos(gen_delta)
    iq[load.bus] -= @. imag(i_load)
    # iq[fault.bus] -= @. fault_iq
    iq[:] += incidence_matrix * line_iq

    # B_mat = build_B_matrix(models)
    # v_dp = @. bus_vd + 1im * bus_vq
    # i_cap = 1im * B_mat * v_dp
    # i_cap_d = zeros(T, length(bus.idx))
    # i_cap_q = zeros(T, length(bus.idx))
    # i_cap_d = real(i_cap)
    # i_cap_q = imag(i_cap)
    # id[:] -= i_cap_d
    # iq[:] -= i_cap_q

    # id[:] += B_mat * bus_vq
    # iq[:] -= B_mat * bus_vd

    omega = 2*pi*60

    w1 = @. omega*C_eq * bus_vq
    w2 = @. omega*C_eq * bus_vd

    id[:] += w1
    iq[:] -= w2

    du[address["balance_d"]] = @. id[non_slack_buses]
    du[address["balance_q"]] = @. iq[non_slack_buses]
end


function solve_dynamic_sim!(du, u, p, t)
    balance!(du, u, p)
end


incidence_matrix = build_incidence_matrix(models)


du = zeros(2*length(non_slack_buses))

begin
    u0 = zeros(length(du))
    u0[address["balance_d"]] = models.bus.vd[non_slack_buses]
    u0[address["balance_q"]] = models.bus.vq[non_slack_buses]
end

mass_matrix = zeros(length(du), length(du))

b_mat = build_B_matrix(models)
C_mat = b_mat ./ (2*pi*60)
C_eq = diag(C_mat)
# add a small value to the diagonal to avoid singularity
# for i in eachindex(C_eq)
#     C_eq[i] = C_eq[i] == 0.0 ? 1e-9 : C_eq[i]
# end

begin

    # vd and vq bus balances
    # i = 1
    # for idx in zip(collect(address["balance_d"]), collect(address["balance_q"]))
    #     mass_matrix[idx[1], idx[1]] = C_eq[i]
    #     mass_matrix[idx[2], idx[2]] = C_eq[i]
    #     i += 1
    # end

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
    
end

p = (address, models, incidence_matrix, C_eq)

tspan = (0.0, 5.0)
prob0 = ODEFunction(solve_dynamic_sim!, mass_matrix=mass_matrix)
prob = ODEProblem(prob0, u0, tspan, p)

# sol = solve(prob, Trapezoid(), adaptive=true], dt = 50e-5, reltol=1e-6, abstol=1e-6)
sol = solve(prob, Rodas5(), adaptive=false, dt = 50e-5)
# sol = solve(prob, Trapezoid(), adaptive=false, dt = 50e-3)

# dt = 50e-5
# prob = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, tspan, p, mass_matrix)
# sol = MyDiffEq.Solve(prob, dt, method=:Euler, adaptive=false);

using Plots

# plot(sol, idxs = address["balance_d"])
# plot(sol, idxs = address["balance_q"])


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
u15 = [sol_u[15] for sol_u in sol.u]
u16 = [sol_u[16] for sol_u in sol.u]
u17 = [sol_u[17] for sol_u in sol.u]
u18 = [sol_u[18] for sol_u in sol.u]
u19 = [sol_u[19] for sol_u in sol.u]
u20 = [sol_u[20] for sol_u in sol.u]
u21 = [sol_u[21] for sol_u in sol.u]
u22 = [sol_u[22] for sol_u in sol.u]
u23 = [sol_u[23] for sol_u in sol.u]
u24 = [sol_u[24] for sol_u in sol.u]
u25 = [sol_u[25] for sol_u in sol.u]
u26 = [sol_u[26] for sol_u in sol.u]

v1 = @. u1 + 1im * u14
v2 = @. u2 + 1im * u15
v3 = @. u3 + 1im * u16
v4 = @. u4 + 1im * u17
v5 = @. u5 + 1im * u18
v6 = @. u6 + 1im * u19
v7 = @. u7 + 1im * u20
v8 = @. u8 + 1im * u21
v9 = @. u9 + 1im * u22
v10 = @. u10 + 1im * u23
v11 = @. u11 + 1im * u24
v12 = @. u12 + 1im * u25
v13 = @. u13 + 1im * u26
plot(abs.(v1))
plot(abs.(v2))
plot(abs.(v3))
plot(abs.(v4))
plot(abs.(v5))
plot(abs.(v6))
plot(abs.(v7))
plot(abs.(v8))
plot(abs.(v9))
plot(abs.(v10))
plot(abs.(v11))
plot(abs.(v12))
plot(abs.(v13))

