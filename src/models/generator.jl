module GeneratorModel

export Generator, solve_generator!

mutable struct Generator{T<:Real}
    idx :: Vector{Int32}
    bus :: Vector{Int32}
    delta :: Vector{T}
    omega :: Vector{T}
    M :: Vector{T}
    p_m :: Vector{T}
    q_m :: Vector{T}
    i_d :: Vector{T}
    i_q :: Vector{T}
    x_d_prime :: Vector{T}
    e_q_prime :: Vector{T}
    d :: Vector{T}

    function Generator{T}(n::Integer) where {T<:Real}
        new{T}(Vector{Int32}(undef, n),
               Vector{Int32}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n))
    end
end

function solve_generator!(du, u, p)
    T = eltype(u)

    Ω = T(2*pi*60)
    d = T(1.0)

    address, models, _, _, non_slack_buses, _ = p

    bus = models.bus
    generator = models.generator

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

end # module
