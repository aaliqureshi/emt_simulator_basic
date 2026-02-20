module BusModel

export Bus, balance!, phasor2DP!

mutable struct Bus{T<:Real}
    idx :: Vector{Int32}
    v :: Vector{T}
    theta :: Vector{T}
    vd:: Vector{T}
    vq:: Vector{T}

    function Bus{T}(n::Integer) where {T<:Real}
        new{T}(Vector{Int32}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n))
    end
end

function phasor2DP!(bus)
    vdq = @. bus.v * exp(1im * bus.theta)
    bus.vd = real(vdq)
    bus.vq = imag(vdq)
end

function balance!(du, u, p)
    T = eltype(u)
    address, models, incidence_matrix, C_eq, non_slack_buses = p

    bus = models.bus
    line = models.line
    generator = models.generator
    fault = models.fault
    load = models.load
    slack = models.slack

    bus_vd = Vector{T}(bus.vd)
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

    # i_load = Y * V  :: constant admittance load model
    # i_load = @. load.y * complex(bus_vd[load.bus], bus_vq[load.bus])
    i_load = @. conj(complex(load.p, load.q) / complex(bus_vd[load.bus], bus_vq[load.bus]))

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

end # module
