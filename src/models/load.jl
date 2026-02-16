module LoadModel

export Load, compute_load_currents

mutable struct Load{T<:Real}
    bus::Vector{Int32}
    p::Vector{T}
    q::Vector{T}

    function Load{T}(n::Integer) where {T<:Real}
        new{T}(Vector{Int32}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n))
    end
end

function compute_load_currents(models)
    v_dp = @. models.bus.vd[models.load.bus] - 1im * models.bus.vq[models.load.bus]
    i_load = @. (models.load.p - 1im * models.load.q) / v_dp
    return real(i_load), imag(i_load)
end

end # module
