module Models

export Bus, Line, Generator, Fault, Load, Shunt, Slack


mutable struct Bus{T<:Real}
    idx :: Vector{Int32}
    v :: Vector{T}
    theta :: Vector{T}
    vd:: Vector{T}
    vq:: Vector{T}
    i_d:: Vector{T}
    i_q:: Vector{T}

    function Bus{T}(n::Integer) where {T<:Real}
        new{T}(Vector{Int32}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n))
    end
end

mutable struct Line{T<:Real}
    idx :: Vector{Int32}
    bus1_idx :: Vector{Int32}
    bus2_idx :: Vector{Int32}
    R :: Vector{T}
    X :: Vector{T}
    i_d:: Vector{T}
    i_q:: Vector{T}

    function Line{T}(n::Integer) where {T<:Real}
        new{T}(Vector{Int32}(undef, n),
               Vector{Int32}(undef, n),
               Vector{Int32}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n), 
               Vector{T}(undef, n))
    end
end

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
               Vector{T}(undef, n))
    end
end

mutable struct Fault{T<:Real}   
    bus::Vector{Int32}
    r_s::Vector{T}
    l_s::Vector{T}

    function Fault{T}(n::Integer) where {T<:Real}
        new{T}(Vector{Int32}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n))
    end
end

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

mutable struct Slack
    bus::Vector{Int32}

    function Slack(n::Integer)
        new(Vector{Int32}(undef, n))
    end
end

end #module