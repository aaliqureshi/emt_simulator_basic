module Models

export Bus, Line, Generator, Fault, Load


mutable struct Bus{T<:Real}
    idx :: Vector{Int32}
    v :: Vector{T}
    theta :: Vector{T}
    vd:: Vector{T}
    vq:: Vector{T}
    i_d:: Vector{T}
    i_q:: Vector{T}
end

mutable struct Line{T<:Real}
    idx :: Vector{Int32}
    bus1_idx :: Vector{Int32}
    bus2_idx :: Vector{Int32}
    R :: Vector{T}
    X :: Vector{T}
    i_d:: Vector{T}
    i_q:: Vector{T}
end

mutable struct Generator{T<:Real}
    idx :: Vector{Int32}
    bus :: Vector{Int32}
    delta :: Vector{T}
    omega :: Vector{T}
    M :: Vector{T}
    p_m :: Vector{T}
    i_d :: Vector{T}
    i_q :: Vector{T}
    x_d_prime :: Vector{T}
    e_q_prime :: Vector{T}
end

mutable struct Fault{T<:Real}
    bus::Int32
    r_s::Float64
    l_s::Float64
    i_d::Vector{T}
    i_q::Vector{T}
end

mutable struct Load{T<:Real}
    bus::Vector{Int32}
    p::Vector{T}
    q::Vector{T}
    i_d::Vector{T}
    i_q::Vector{T}
end

end