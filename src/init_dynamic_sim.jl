# TODO: add load currents in the equations


n_gens = length(models.generator.bus)
n_lines = length(models.line.idx)

idx_delta = 1:n_gens
idx_omega = idx_delta[end]+1 : idx_delta[end]+n_gens
idx_e_q_prime = idx_omega[end]+1 : idx_omega[end]+n_gens
idx_gen_id = idx_e_q_prime[end]+1 : idx_e_q_prime[end]+n_gens
idx_gen_iq = idx_gen_id[end]+1 : idx_gen_id[end]+n_gens
idx_line_id = idx_gen_iq[end]+1 : idx_gen_iq[end]+n_lines
idx_line_iq = idx_line_id[end]+1 : idx_line_id[end]+n_lines

address = Dict(
    "delta" => idx_delta,
    "omega" => idx_omega,
    "e_q_prime" => idx_e_q_prime,
    "gen_id" => idx_gen_id,
    "gen_iq" => idx_gen_iq,
    "line_id" => idx_line_id,
    "line_iq" => idx_line_iq
)


function init_line!(du, u, address, models)
    T = eltype(u)

    bus = models.bus
    line = models.line
    generator = models.generator
    load = models.load
    slack = models.slack

    gen_delta = u[address["delta"]]
    gen_omega = u[address["omega"]]
    gen_e_q_prime = u[address["e_q_prime"]]
    gen_id = u[address["gen_id"]]
    gen_iq = u[address["gen_iq"]]
    line_id = u[address["line_id"]]
    line_iq = u[address["line_iq"]]

    bus_vd = Vector{T}(bus.vd)  # converts Float64 -> T safely
    bus_vq = Vector{T}(bus.vq)

    du[address["line_id"]] = @. ((bus_vd[line.bus1_idx] - bus_vd[line.bus2_idx] - line.R * line_id) + (line.X * line_iq))
    du[address["line_iq"]] = @. ((bus_vq[line.bus1_idx] - bus_vq[line.bus2_idx] - line.R * line_iq) - (line.X * line_id))

end

function init_gens!(du, u, address, models)
    T = eltype(u)

    Ω = T(2*pi*60)
    d = T(1.0)

    bus = models.bus
    generator = models.generator
    line = models.line
    load = models.load
    slack = models.slack

    gen_delta = u[address["delta"]]
    gen_omega = u[address["omega"]]
    gen_e_q_prime = u[address["e_q_prime"]]
    gen_id = u[address["gen_id"]]
    gen_iq = u[address["gen_iq"]]
    line_id = u[address["line_id"]]
    line_iq = u[address["line_iq"]]

    bus_vd = Vector{T}(bus.vd)  # converts Float64 -> T safely
    bus_vq = Vector{T}(bus.vq)

    du[address["delta"]] = @. Ω * (gen_omega - one(T))
    du[address["omega"]] = @. ((generator.p_m - (gen_id * bus_vd[generator.bus] * sin(gen_delta) - 
                                     gen_id * bus_vq[generator.bus] * cos(gen_delta) + 
                                     gen_iq * bus_vd[generator.bus] * cos(gen_delta) + 
                                     gen_iq * bus_vq[generator.bus] * sin(gen_delta)) - 
                                     d * (gen_omega - one(T))) / (generator.M))
    du[address["gen_id"]] = @. gen_e_q_prime - gen_id * generator.x_d_prime - 
                                bus_vd[generator.bus] * cos(gen_delta) - 
                                bus_vq[generator.bus] * sin(gen_delta)
    du[address["gen_iq"]] = @. gen_iq * generator.x_d_prime - 
                                bus_vd[generator.bus] * sin(gen_delta) + 
                                bus_vq[generator.bus] * cos(gen_delta)

    ## temporary, remove later
    # du[address["e_q_prime"]] = @. gen_id - line_iq[1]
    du[address["e_q_prime"]] = @. line_id[1] - (
                                  gen_id * cos(gen_delta - pi/2) - 
                                  gen_iq * sin(gen_delta - pi/2))

end

function init_system!(u, p)
    T = eltype(u)
    du = zeros(T, length(u))
    address, models = p

    init_gens!(du, u, address, models)
    init_line!(du, u, address, models)

    return du
end

p = (address, models)
du = zeros(5 * n_gens + 2 * n_lines)
u0 = ones(length(du))

prob = NonlinearProblem(init_system!, u0, p)
sol = solve(prob, abstol = 1e-8, reltol = 1e-8)



