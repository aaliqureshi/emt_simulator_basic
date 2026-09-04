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
    # r_eff = (1e10)^(1-lambda) * (fault.r_fault[1])^(lambda)
    # r_eff = (1e10)^(lambda) * (fault.r_fault[1])^(1-lambda) 
    # g_fault = 1 / r_eff

    # du[address["fault_id"]] = @. g_fault * bus_vd[fault.bus] - fault_id
    # du[address["fault_iq"]] = @. g_fault * bus_vq[fault.bus] - fault_iq

    # reactive fault
    x_eff = (1e10)^(1-lambda) * (fault.x_fault[1])^(lambda)
    # x_eff = (1e10)^(lambda) * (fault.x_fault[1])^(1-lambda) 
    # x_eff = t < 2.0 ? (1e10)^(1-lambda) * (fault.x_fault[1])^(lambda) : 1e10
    # x_eff = 1.0 <= t <= 1.1 ? fault.x_fault[1] : 1e10
    # x_eff = models.fault.x_fault[1]
    b_fault = 1 / x_eff
    du[address["fault_id"]] = @. b_fault * bus_vd[fault.bus] + fault_iq
    du[address["fault_iq"]] = @. b_fault * bus_vq[fault.bus] - fault_id

    # resistive + reactive fault
    # r_eff = (1e10)^(lambda) * (fault.r_fault[1])^(1-lambda) 
    # x_eff = (1e10)^(lambda) * (fault.x_fault[1])^(1-lambda) 
    # r_eff = (1e10)^(1-lambda) * (fault.r_fault[1])^(lambda) 
    # x_eff = (1e10)^(1-lambda) * (fault.x_fault[1])^(lambda) 
    # g_fault = 1 / r_eff
    # b_fault = 1 / x_eff
    # du[address["fault_id"]] = @. g_fault * b_fault * bus_vd[fault.bus] - b_fault * fault_id + g_fault * fault_iq
    # du[address["fault_iq"]] = @. g_fault * b_fault * bus_vq[fault.bus] - g_fault * fault_id - b_fault * fault_iq
end

end # module
