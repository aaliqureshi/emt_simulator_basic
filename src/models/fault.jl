module FaultModel

export Fault, solve_fault!

mutable struct Fault{T<:Real}
    bus::Vector{Int32}
    r_fault::Vector{T}
    r_open::Vector{T}
    x_fault::Vector{T}
    x_open::Vector{T}
    l_fault::Vector{T}
    l_open::Vector{T}
    t_fault::Vector{T}
    t_clear::Vector{T}

    function Fault{T}(n::Integer) where {T<:Real}
        new{T}(Vector{Int32}(undef, n),
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

function solve_fault!(du, u, p, t)
    T = eltype(u)
    Ω = T(2*pi*60)

    address, models, _, _, non_slack_buses, lambda = p

    bus = models.bus
    fault = models.fault

    fault_id = u[address["fault_id"]]
    fault_iq = u[address["fault_iq"]]

    bus_vd = Vector{T}(bus.vd)
    bus_vq = Vector{T}(bus.vq)

    bus_vd[non_slack_buses] = @. u[address["balance_d"]]
    bus_vq[non_slack_buses] = @. u[address["balance_q"]]

    # r_eff = fault.t_fault[1] <= t <= fault.t_clear[1] ? fault.r_fault[1] : fault.r_open[1]
    r_eff = (1e10)^(1-lambda) * (fault.r_fault[1])^(lambda)
    g_fault = 1 / r_eff
    #TODO: add reactive fault 

    du[address["fault_id"]] = @. g_fault * bus_vd[fault.bus] - fault_id
    du[address["fault_iq"]] = @. g_fault * bus_vq[fault.bus] - fault_iq
end

end # module
