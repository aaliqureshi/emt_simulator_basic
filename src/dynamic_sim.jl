module DynamicSim

export solve_dynamic_sim!, build_dynamic_address, build_mass_matrix, build_initial_conditions

include("Models.jl"); using .Models
# include("algebraic_solvers.jl")

function solve_dynamic_sim!(du, u, p, t)
    solve_generator!(du, u, p)
    solve_line!(du, u, p, t)
    solve_fault!(du, u, p, t)
    balance!(du, u, p)
end

function build_dynamic_address(sys)
    models = sys.models
    non_slack_buses = sys.non_slack_buses

    n_gens = length(models.generator.bus)
    n_lines = length(models.line.idx)
    n_buses = length(models.bus.idx)
    n_faults = length(models.fault.bus)

    idx_delta = 1:n_gens
    idx_omega = idx_delta[end]+1 : idx_delta[end]+n_gens
    idx_line_id = idx_omega[end]+1 : idx_omega[end]+n_lines
    idx_line_iq = idx_line_id[end]+1 : idx_line_id[end]+n_lines
    idx_balance_d = idx_line_iq[end]+1 : idx_line_iq[end]+length(non_slack_buses)
    idx_balance_q = idx_balance_d[end]+1 : idx_balance_d[end]+length(non_slack_buses)
    idx_fault_id = idx_balance_q[end]+1 : idx_balance_q[end]+n_faults
    idx_fault_iq = idx_fault_id[end]+1 : idx_fault_id[end]+n_faults
    idx_gen_id = idx_fault_iq[end]+1 : idx_fault_iq[end]+n_gens
    idx_gen_iq = idx_gen_id[end]+1 : idx_gen_id[end]+n_gens

    address = Dict(
        "delta" => idx_delta,
        "omega" => idx_omega,
        "gen_id" => idx_gen_id,
        "gen_iq" => idx_gen_iq,
        "line_id" => idx_line_id,
        "line_iq" => idx_line_iq,
        "fault_id" => idx_fault_id,
        "fault_iq" => idx_fault_iq,
        "balance_d" => idx_balance_d,
        "balance_q" => idx_balance_q
    )

    return address
end

function build_mass_matrix(sys, address)
    models = sys.models
    non_slack_buses = sys.non_slack_buses
    C_eq = sys.C_eq

    n_gens = length(models.generator.bus)
    n_lines = length(models.line.idx)
    n_faults = length(models.fault.bus)
    n_states = 4*n_gens + 2*n_lines + 2*n_faults + 2*length(non_slack_buses)

    mass_matrix = zeros(n_states, n_states)

    # bus balance equations (capacitance)
    for (i, bus_idx) in enumerate(address["balance_d"])
        real_bus_idx = non_slack_buses[i]
        mass_matrix[bus_idx, bus_idx] = C_eq[real_bus_idx]
        q_idx = address["balance_q"][i]
        mass_matrix[q_idx, q_idx] = C_eq[real_bus_idx]
    end

    # line id (inductance)
    i = 1
    for idx in collect(address["line_id"])
        mass_matrix[idx, idx] = models.line.L[i]
        i += 1
    end

    # line iq (inductance)
    i = 1
    for idx in collect(address["line_iq"])
        mass_matrix[idx, idx] = models.line.L[i]
        i += 1
    end

    # # delta
    for i in collect(address["delta"])
        mass_matrix[i,i] = 1.0
    end

    # omega (inertia)
    idx = 1
    for i in collect(address["omega"])
        mass_matrix[i,i] = models.generator.M[idx]
        idx += 1
    end

    return mass_matrix
end

function build_initial_conditions(sys, address)
    models = sys.models
    non_slack_buses = sys.non_slack_buses

    n_gens = length(models.generator.bus)
    n_lines = length(models.line.idx)
    n_faults = length(models.fault.bus)
    n_states = 4*n_gens + 2*n_lines + 2*n_faults + 2*length(non_slack_buses)

    u0 = zeros(n_states)
    u0[address["delta"]] = models.generator.delta
    u0[address["omega"]] = models.generator.omega
    u0[address["gen_id"]] = models.generator.i_d
    u0[address["gen_iq"]] = models.generator.i_q
    u0[address["line_id"]] = models.line.i_d
    u0[address["line_iq"]] = models.line.i_q
    u0[address["fault_id"]] = [0.0]
    u0[address["fault_iq"]] = [0.0]
    u0[address["balance_d"]] = models.bus.vd[non_slack_buses]
    u0[address["balance_q"]] = models.bus.vq[non_slack_buses]

    return u0
end

end # module
