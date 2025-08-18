module PowerFlow

export PowerFlowSystem, solve_powerflow, phasor2DP!, run_powerflow_example

using Pkg
Pkg.activate(".")

include("Models.jl")
using .Models

include("data_loader.jl")
using .DataLoader

using OrdinaryDiffEq

"""
    PowerFlowSystem

A struct to hold all the components and state information for a power flow simulation.
"""
struct PowerFlowSystem{T<:Real}
    models::NamedTuple
    address::Dict{String, UnitRange{Int}}
    u0::Vector{T}
    tspan::Tuple{Float64, Float64}
    vd_update_buses::Vector{Int}
    vq_update_buses::Vector{Int}
    non_slack_buses::Vector{Int}
    
    function PowerFlowSystem(data_file::String)
        # Load data
        models_dict = load_data(data_file)
        models = (bus=models_dict["bus"], 
                  line=models_dict["line"], 
                  generator=models_dict["generator"], 
                  fault=models_dict["fault"], 
                  load=models_dict["load"], 
                  slack=models_dict["slack"])
        
        # Initialize bus voltages
        phasor2DP!(models.bus)
        
        # Convert line reactance to per unit
        models.line.X = models.line.X ./ (2*pi*60)
        
        # Calculate dimensions
        n_gen = length(models.generator.bus)
        n_line = length(models.line.idx)
        n_bus = length(models.bus.idx)
        n_load = length(models.load.bus)
        n_slack = length(models.slack.bus)
        
        # Define update buses
        vd_update_buses = sort(setdiff(models.load.bus, models.generator.bus))
        vq_update_buses = sort(union(models.load.bus, models.generator.bus))
        non_slack_buses = setdiff(models.bus.idx, models.slack.bus)
        
        # Create state indexing
        idx_line_d = 1:n_line
        idx_line_q = idx_line_d[end]+1 : idx_line_d[end]+n_line
        idx_balance_q = idx_line_q[end]+1 : idx_line_q[end]+(n_load)
        idx_balance_d = idx_balance_q[end]+1 : idx_balance_q[end]+(n_load + n_gen)
        
        address = Dict(
            "line_d" => idx_line_d,
            "line_q" => idx_line_q,
            "balance_d" => idx_balance_d,
            "balance_q" => idx_balance_q
        )
        
        # Initial state vector
        u0 = vcat(models.line.i_d, models.line.i_q, 
                  models.bus.v[vd_update_buses], models.bus.theta[vq_update_buses])
        
        tspan = (0.0, 10.0)
        
        new{eltype(u0)}(models, address, u0, tspan, vd_update_buses, vq_update_buses, non_slack_buses)
    end
end

"""
    phasor2DP!(bus)

Convert phasor representation (magnitude and angle) to d-q components.
"""
function phasor2DP!(bus)
    vdq = @. bus.v * exp(1im * bus.theta)
    bus.vd = real(vdq)
    bus.vq = imag(vdq)
end

"""
    line_equation!(du, u, address, models, vd_update_buses, vq_update_buses)

Compute the line differential equations.
"""
function line_equation!(du, u, address::Dict, models::NamedTuple, 
                       vd_update_buses::Vector{Int}, vq_update_buses::Vector{Int})
    T = eltype(u)
    Ω = T(2*pi*60)

    bus = models.bus
    line = models.line

    line_id = u[address["line_d"]]
    line_iq = u[address["line_q"]]

    # Update bus voltages from state vector
    bus_v = Vector{T}(bus.v)
    bus_theta = Vector{T}(bus.theta)
    bus_v[vd_update_buses] .= u[address["balance_q"]]
    bus_theta[vq_update_buses] .= u[address["balance_d"]]

    # Convert to d-q components
    bus_vd = @. bus_v * cos(bus_theta)
    bus_vq = @. bus_v * sin(bus_theta)

    # Line equations
    du[address["line_d"]] = @. ((bus_vd[line.bus1_idx] - bus_vd[line.bus2_idx] - line.R * line_id) / line.X) + (Ω * line_iq)
    du[address["line_q"]] = @. ((bus_vq[line.bus1_idx] - bus_vq[line.bus2_idx] - line.R * line_iq) / line.X) - (Ω * line_id)
end

"""
    power_balance!(du, u, address, models, vd_update_buses, vq_update_buses, non_slack_buses)

Compute the power balance equations.
"""
function power_balance!(du, u, address::Dict, models::NamedTuple, 
                       vd_update_buses::Vector{Int}, vq_update_buses::Vector{Int}, 
                       non_slack_buses::Vector{Int})
    T = eltype(u)

    bus = models.bus
    line = models.line
    generator = models.generator
    load = models.load

    line_id = u[address["line_d"]]
    line_iq = u[address["line_q"]]

    # Update bus voltages from state vector
    bus_v = Vector{T}(bus.v)
    bus_theta = Vector{T}(bus.theta)
    bus_v[vd_update_buses] .= u[address["balance_q"]]
    bus_theta[vq_update_buses] .= u[address["balance_d"]]

    # Convert to d-q components
    bus_vd = @. bus_v * cos(bus_theta)
    bus_vq = @. bus_v * sin(bus_theta)

    # Initialize power arrays
    real_power = zeros(T, length(bus.idx))
    reactive_power = zeros(T, length(bus.idx))

    # Load and generator powers
    p_load = load.p
    q_load = load.q
    p_gen = generator.p_m
    q_gen = generator.q_m

    # Line power flows
    p_line_from = @. bus_vd[line.bus1_idx] * line_id + bus_vq[line.bus1_idx] * line_iq
    q_line_from = @. bus_vq[line.bus1_idx] * line_id - bus_vd[line.bus1_idx] * line_iq
    p_line_to = @. bus_vd[line.bus2_idx] * line_id + bus_vq[line.bus2_idx] * line_iq
    q_line_to = @. bus_vq[line.bus2_idx] * line_id - bus_vd[line.bus2_idx] * line_iq

    # Power balance equations
    real_power[load.bus] -= p_load
    real_power[generator.bus] += p_gen
    real_power[line.bus1_idx] -= p_line_from
    real_power[line.bus2_idx] += p_line_to

    reactive_power[load.bus] -= q_load
    reactive_power[generator.bus] += q_gen
    reactive_power[line.bus1_idx] -= q_line_from
    reactive_power[line.bus2_idx] += q_line_to

    # Set derivatives
    du[address["balance_d"]] = @. real_power[non_slack_buses]
    du[address["balance_q"]] = @. reactive_power[load.bus]
end

"""
    powerflow!(du, u, p, t)

Main power flow differential equation system.
"""
function powerflow!(du, u, p, t)
    address, models, vd_update_buses, vq_update_buses, non_slack_buses = p
    
    # Compute line equations
    line_equation!(du, u, address, models, vd_update_buses, vq_update_buses)
    
    # Compute power balance equations
    power_balance!(du, u, address, models, vd_update_buses, vq_update_buses, non_slack_buses)
end

"""
    solve_powerflow(system::PowerFlowSystem; dt=50e-5)

Solve the power flow system and return the solution.
"""
function solve_powerflow(system::PowerFlowSystem; dt=50e-5)
    # Create mass matrix (identity for this case)
    M = zeros(length(system.u0), length(system.u0))
    
    # Create problem parameters
    p = (system.address, system.models, system.vd_update_buses, 
         system.vq_update_buses, system.non_slack_buses)
    
    # Create and solve ODE problem
    prob0 = ODEFunction(powerflow!, mass_matrix=M)
    prob = ODEProblem(prob0, system.u0, system.tspan, p)
    sol = solve(prob, Trapezoid(), adaptive=false, dt=dt, reltol=1e-10, abstol=1e-10)
    
    return sol
end

"""
    update_models!(system::PowerFlowSystem, sol)

Update the model states with the final solution values.
"""
function update_models!(system::PowerFlowSystem, sol)
    final_state = sol.u[end]
    
    # Update line currents
    system.models.line.i_d .= final_state[system.address["line_d"]]
    system.models.line.i_q .= final_state[system.address["line_q"]]
    
    # Update bus voltages and angles
    system.models.bus.v[system.vd_update_buses] .= final_state[system.address["balance_q"]]
    system.models.bus.theta[system.vq_update_buses] .= final_state[system.address["balance_d"]]
    
    # Update d-q components
    phasor2DP!(system.models.bus)
end

# Example usage function
"""
    run_powerflow_example(data_file="cases/SMIB_Chow/SMIB_RL_Line_DrCui.xlsx")

Run a complete power flow simulation example.
"""
function run_powerflow_example(data_file="cases/SMIB_Chow/SMIB_RL_Line_DrCui.xlsx")
    # Create power flow system
    system = PowerFlowSystem(data_file)
    
    # Solve the system
    sol = solve_powerflow(system)
    
    # Update models with final solution
    update_models!(system, sol)
    
    # Print results
    println("Power Flow Solution Complete")
    println("Line currents (d-axis): ", system.models.line.i_d)
    println("Line currents (q-axis): ", system.models.line.i_q)
    println("Bus voltages: ", system.models.bus.v)
    println("Bus angles: ", system.models.bus.theta)
    
    return system, sol
end

end # module