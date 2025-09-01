# run after powerflow_traditional.jl
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


function init_line!(du, u, p)
    T = eltype(u)

    address, models, _ = p

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

function init_gens!(du, u, p)
    T = eltype(u)

    Ω = T(2*pi*60)
    d = T(1.0)

    address, models, incidence_matrix = p

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
    # TODO: add full KCL equation here: 
    ## gen_id - load_id + network_id = 0
    du[address["e_q_prime"]] = @. line_id[1] - (
                                  gen_id * cos(gen_delta - pi/2) - 
                                  gen_iq * sin(gen_delta - pi/2))

end

function init_system!(u, p)
    T = eltype(u)
    du = zeros(T, length(u))
    # address, models = p

    init_gens!(du, u, p)
    init_line!(du, u, p)

    return du
end

## add incidence matrix here
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


incidence_matrix = create_incidence_matrix(models)


p = (address, models, incidence_matrix)
du = zeros(5 * n_gens + 2 * n_lines)
u0 = ones(length(du))

prob = NonlinearProblem(init_system!, u0, p)
sol = solve(prob, abstol = 1e-8, reltol = 1e-8)


## update models

models.generator.delta .= sol.u[address["delta"]]
models.generator.omega .= sol.u[address["omega"]]
models.generator.e_q_prime .= sol.u[address["e_q_prime"]]
models.generator.i_d .= sol.u[address["gen_id"]]
models.generator.i_q .= sol.u[address["gen_iq"]]
models.line.i_d .= sol.u[address["line_id"]]
models.line.i_q .= sol.u[address["line_iq"]]


