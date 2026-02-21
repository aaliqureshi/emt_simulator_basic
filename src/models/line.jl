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

function solve_line!(du, u, p, t)
    T = eltype(u)
    Ω = T(2*pi*60)

    address, models, _, _, non_slack_buses, lambda = p

    bus = models.bus
    line = models.line

    line_id = u[address["line_id"]]
    line_iq = u[address["line_iq"]]

    bus_vd = Vector{T}(bus.vd)
    bus_vq = Vector{T}(bus.vq)

    bus_vd[non_slack_buses] = @. u[address["balance_d"]]
    bus_vq[non_slack_buses] = @. u[address["balance_q"]]

    # models.line.R[2] = 1.0 <= t <= 1.1 ? 1e10 : 0.05403
    # models.line.X[2] = 1.0 <= t <= 1.1 ? 1e10 : 0.22304

    # models.line.R[2] = t > 2.0 ? 1e10 : 0.05403
    # models.line.X[2] = t > 2.0 ? 1e10 : 0.22304

    # lines = [2, 5, 7]
    lines = [2]
    # models.line.R[lines] .= t > 2.0 ? 1e10 : models.line.R[lines]
    # models.line.X[lines] .= t > 2.0 ? 1e10 : models.line.X[lines]

    # models.line.R[lines] .= 1e10
    # models.line.X[lines] .= 1e10

    # models.line.R[2] = (0.05403)^(1-lambda) * (1e10)^(lambda)
    # models.line.X[2] = (0.22304)^(1-lambda) * (1e10)^(lambda)


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
