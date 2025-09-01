# run after init_dynamic_sim.jl

using OrdinaryDiffEq


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
idx_gen_id = idx_fault_iq[end]+1 : idx_fault_iq[end]+n_gens
idx_gen_iq = idx_gen_id[end]+1 : idx_gen_id[end]+n_gens
idx_balance_d = idx_gen_iq[end]+1 : idx_gen_iq[end]+(length(non_slack_buses))
idx_balance_q = idx_balance_d[end]+1 : idx_balance_d[end]+(length(non_slack_buses))




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


## create incidence matrix

function create_incidence_matrix(models)
    bus1_idx = models.line.bus1_idx
    bus2_idx = models.line.bus2_idx

    num_bus = length(models.bus.idx)
    num_line = length(models.line.idx)
    incidence_matrix = zeros(Int32, num_bus, num_line)

    for line in collect(1:num_line)
        incidence_matrix[bus1_idx[line], line] = -1
        incidence_matrix[bus2_idx[line], line] = 1
    end

    return incidence_matrix
end

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
    du[address["omega"]] = @. (generator.p_m - (gen_id * bus_vd[generator.bus] * sin(gen_delta) - 
                                     gen_id * bus_vq[generator.bus] * cos(gen_delta) + 
                                     gen_iq * bus_vd[generator.bus] * cos(gen_delta) + 
                                     gen_iq * bus_vq[generator.bus] * sin(gen_delta)) - 
                                     d * (gen_omega - one(T))) / (generator.M)
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

    L = @. T(line.X) / Ω
    

    bus_vd = Vector{T}(bus.vd)  # converts Float64 -> T safely
    bus_vq = Vector{T}(bus.vq)

    bus_vd[non_slack_buses] .= u[address["balance_d"]]
    bus_vq[non_slack_buses] .= u[address["balance_q"]]

    du[address["line_id"]] = @. ((bus_vd[line.bus1_idx] - bus_vd[line.bus2_idx] - line.R * line_id) / L + (Ω * line_iq))
    du[address["line_iq"]] = @. ((bus_vq[line.bus1_idx] - bus_vq[line.bus2_idx] - line.R * line_iq) / L - (Ω * line_id))

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

    Ls = @. T(fault.l_s * 1e100) / Ω

    bus_vd = Vector{T}(bus.vd)  # converts Float64 -> T safely
    bus_vq = Vector{T}(bus.vq)

    bus_vd[non_slack_buses] .= u[address["balance_d"]]
    bus_vq[non_slack_buses] .= u[address["balance_q"]]

    du[address["fault_id"]] = @. ((bus_vd[fault.bus] - fault.r_s * fault_id) / Ls + (Ω * fault_iq))
    du[address["fault_iq"]] = @. ((bus_vq[fault.bus] - fault.r_s * fault_iq) / Ls - (Ω * fault_id))

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

    i_load = @. (load.p - 1im * load.q) / (bus_vd[load.bus] + 1im * bus_vq[load.bus])

    id[generator.bus] += @. gen_id * cos(gen_delta - pi/2) - gen_iq * sin(gen_delta - pi/2)
    id[load.bus] -= @. real(i_load)
    id[fault.bus] -= @. fault_id
    id[:] += incidence_matrix * line_id

    iq[generator.bus] += @. gen_id * sin(gen_delta - pi/2) + gen_iq * cos(gen_delta - pi/2)
    iq[load.bus] -= @. imag(i_load)
    iq[fault.bus] -= @. fault_iq
    iq[:] += incidence_matrix * line_iq

    du[address["balance_d"]] = @. id[non_slack_buses]
    du[address["balance_q"]] = @. iq[non_slack_buses]
end


function solve_dynamic_sim!(du, u, p, t)
    solve_generator!(du, u, p)
    solve_line!(du, u, p)
    solve_fault!(du, u, p)
    balance!(du, u, p)
end


incidence_matrix = create_incidence_matrix(models)
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
    u0[address["fault_id"]] = models.fault.i_d
    u0[address["fault_iq"]] = models.fault.i_q
    u0[address["balance_d"]] = models.bus.vd[non_slack_buses]
    u0[address["balance_q"]] = models.bus.vq[non_slack_buses]
end

mass_matrix = zeros(length(du), length(du))
for i in range(1, (2*n_gens + 2*n_lines + 2*n_faults))
    mass_matrix[i,i] = 1.0
end

tspan = (0.0, 10.0)
prob0 = ODEFunction(solve_dynamic_sim!, mass_matrix=mass_matrix)
prob = ODEProblem(prob0, u0, tspan, p)
sol = solve(prob, Trapezoid(), adaptive=false, dt = 50e-5)

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




### quick test

incidence_matrix = create_incidence_matrix(models)
du = zeros(4*n_gens + 2*n_lines + 2*n_faults + 2*length(non_slack_buses))

u0 = vcat(models.generator.delta, models.generator.omega,
          models.line.i_d, models.line.i_q,
          models.fault.i_d, models.fault.i_q,
          models.generator.i_d, models.generator.i_q,
          models.bus.vd[non_slack_buses], models.bus.vq[non_slack_buses])

bus = models.bus
generator = models.generator
line = models.line
load = models.load
slack = models.slack
fault = models.fault

gen_delta = u0[address["delta"]]
gen_omega = u0[address["omega"]]
gen_id = u0[address["gen_id"]]
gen_iq = u0[address["gen_iq"]]
line_id = u0[address["line_id"]]
line_iq = u0[address["line_iq"]]
fault_id = u0[address["fault_id"]]
fault_iq = u0[address["fault_iq"]]
balance_d = u0[address["balance_d"]]
balance_q = u0[address["balance_q"]]

bus_vd = Vector{Float64}(bus.vd)  # converts Float64 -> T safely
bus_vq = Vector{Float64}(bus.vq)

bus_vd[non_slack_buses] .= balance_d
bus_vq[non_slack_buses] .= balance_q

T = Float64
Ω = T(2*pi*60)
d = T(1.0)
L = (line.X)/Ω
Ls = 1e100 .* fault.l_s / Ω

id = zeros(T, length(bus.idx))
iq = zeros(T, length(bus.idx))

i_load = (load.p - 1im * load.q) / (bus_vd[load.bus] + 1im * bus_vq[load.bus])


du[address["delta"]] = @. Ω * (gen_omega - one(T))
du[address["omega"]] = @. (generator.p_m - (gen_id * bus_vd[generator.bus] * sin(gen_delta) - 
                                 gen_id * bus_vq[generator.bus] * cos(gen_delta) + 
                                 gen_iq * bus_vd[generator.bus] * cos(gen_delta) + 
                                 gen_iq * bus_vq[generator.bus] * sin(gen_delta)) - 
                                 d * (gen_omega - one(T))) / (generator.M)
du[address["gen_id"]] = @. generator.e_q_prime - gen_id * generator.x_d_prime - 
                            bus_vd[generator.bus] * cos(gen_delta) - 
                            bus_vq[generator.bus] * sin(gen_delta)
du[address["gen_iq"]] = @. gen_iq * generator.x_d_prime - 
                            bus_vd[generator.bus] * sin(gen_delta) + 
                            bus_vq[generator.bus] * cos(gen_delta)

du[address["line_id"]] = @. ((bus_vd[line.bus1_idx] - bus_vd[line.bus2_idx] - line.R * line_id) / L + (Ω * line_iq))
du[address["line_iq"]] = @. ((bus_vq[line.bus1_idx] - bus_vq[line.bus2_idx] - line.R * line_iq) / L - (Ω * line_id))

du[address["fault_id"]] = @. ((bus_vd[fault.bus] - fault.r_s * fault_id) / Ls + (Ω * fault_iq))
du[address["fault_iq"]] = @. ((bus_vq[fault.bus] - fault.r_s * fault_iq) / Ls - (Ω * fault_id))

id[generator.bus] += @. gen_id * cos(gen_delta - pi/2) - gen_iq * sin(gen_delta - pi/2)
id[load.bus] -= @. real(i_load)
id[fault.bus] -= @. fault_id
id[:] += incidence_matrix * line_id

iq[generator.bus] += @. gen_id * sin(gen_delta - pi/2) + gen_iq * cos(gen_delta - pi/2)
iq[load.bus] -= @. imag(i_load)
iq[fault.bus] -= @. fault_iq
iq[:] += incidence_matrix * line_iq

du[address["balance_d"]] = @. id[non_slack_buses]
du[address["balance_q"]] = @. iq[non_slack_buses]




prob = NonlinearProblem(solve_dynamic_sim!, u0, p)
sol = solve(prob, reltol = 1e-10)


system.models.generator.delta .= sol.u[address["delta"]]