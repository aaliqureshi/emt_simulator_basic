module LineModel

export Line, solve_line!, compute_line_currents!

mutable struct Line{T<:Real}
    idx :: Vector{Int32}
    bus1_idx :: Vector{Int32}
    bus2_idx :: Vector{Int32}
    R :: Vector{T}
    X :: Vector{T}
    B :: Vector{T}
    C :: Vector{T}
    L :: Vector{T}
    i_d:: Vector{T}
    i_q:: Vector{T}

    function Line{T}(n::Integer) where {T<:Real}
        new{T}(Vector{Int32}(undef, n),
               Vector{Int32}(undef, n),
               Vector{Int32}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n))
    end
end

function solve_line!(du, u, p)
    T = eltype(u)
    Ω = T(2*pi*60)

    address, models, _, _, non_slack_buses = p

    bus = models.bus
    line = models.line

    line_id = u[address["line_id"]]
    line_iq = u[address["line_iq"]]

    bus_vd = Vector{T}(bus.vd)
    bus_vq = Vector{T}(bus.vq)

    bus_vd[non_slack_buses] = @. u[address["balance_d"]]
    bus_vq[non_slack_buses] = @. u[address["balance_q"]]

    du[address["line_id"]] = @. bus_vd[line.bus1_idx] - bus_vd[line.bus2_idx] - line.R * line_id + line.X * line_iq
    du[address["line_iq"]] = @. bus_vq[line.bus1_idx] - bus_vq[line.bus2_idx] - line.R * line_iq - line.X * line_id
end

function compute_line_currents!(models)
    v_dp = @. models.bus.vd + 1im * models.bus.vq
    i_line = @. (v_dp[models.line.bus1_idx] - v_dp[models.line.bus2_idx]) / (models.line.R + 1im * models.line.X)
    models.line.i_d .= real(i_line)
    models.line.i_q .= imag(i_line)
end

end # module
