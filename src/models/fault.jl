module FaultModel

export Fault, solve_fault!

mutable struct Fault{T<:Real}
    bus::Vector{Int32}
    r_s::Vector{T}
    x_s::Vector{T}
    l_s::Vector{T}
    tf::Vector{T}
    tc::Vector{T}

    function Fault{T}(n::Integer) where {T<:Real}
        new{T}(Vector{Int32}(undef, n),
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

    address, models, _, _, non_slack_buses = p

    bus = models.bus
    fault = models.fault

    fault_id = u[address["fault_id"]]
    fault_iq = u[address["fault_iq"]]

    Ls = @. T(fault.l_s) / Ω

    bus_vd = Vector{T}(bus.vd)
    bus_vq = Vector{T}(bus.vq)

    bus_vd[non_slack_buses] = @. u[address["balance_d"]]
    bus_vq[non_slack_buses] = @. u[address["balance_q"]]

    r_s = fault.r_s
    r_eff = (1e-2)^(t / 0.1) * (1e10)^(1-t/0.1)

    # g_fault = @. 1 / (r_eff)
    g_fault = 1 / r_s

    du[address["fault_id"]] = @. g_fault * bus_vd[fault.bus] - fault_id
    du[address["fault_iq"]] = @. g_fault * bus_vq[fault.bus] - fault_iq
end

end # module
