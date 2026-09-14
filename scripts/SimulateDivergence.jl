module Simulate
using Barq
using MyDiffEq

export run_simulation

function _get_data_file(case_bus)
    file_dict = Dict(39 => "cases/Fault_Cases/ieee39_fault.xlsx",
                     1 => "cases/Fault_Cases/SMIB_RL_Line_DrCui.xlsx",
                     118 => "cases/Fault_Cases/case118_gc.xlsx",
                     14 => "cases/Fault_Cases/ieee14_fault_barq_no_shunt.xlsx",
                )

    !haskey(file_dict, case_bus) && error("Wrong case file")
    return file_dict[case_bus]
end


function _create_system(data_file, load_factor)
    models = load_data(data_file)
    sys = build_system(models)
    sys.models.load.p[:] .*= load_factor
    return sys
end

function _add_fault!(sys, case_bus)
    if case_bus == 39
        # sys.models.fault.bus = [20]
        # sys.models.fault.x_fault[1] = 0.015
        sys.models.fault.bus = [24]
        # sys.models.fault.x_fault[1] = 0.001
        sys.models.fault.x_fault[1] = 0.006
    elseif case_bus == 118
        sys.models.fault.bus = [12]
        # sys.models.fault.x_fault[1] = 0.01
        sys.models.fault.x_fault[1] = 0.009
        # sys.models.fault.x_fault[1] = 0.05
    end
    nothing
end

function _create_dynamic_data(sys)
    address = build_dynamic_address(sys)
    mass_matrix = build_mass_matrix(sys, address)
    u0 = build_initial_conditions(sys, address)
    p = (address, sys.models, sys.incidence_matrix, sys.C_eq, 
       sys.non_slack_buses, 0.0)
    return u0, p, mass_matrix
end

function _run_sim(u0, p, mass_matrix;
                  dt=5e-4, 
                  lambda=0.0,
                  method=:Euler,
                  always_new=true,
                  tstops=[],
                  adaptive=false,
                  t_end=0.1,
    )
    prob = MyDiffEq.ODEProblem(
                    solve_dynamic_sim!, 
                    u0, 
                    (0.0, t_end), 
                    p, 
                    mass_matrix,
    )
    sol = MyDiffEq.Solve(
                   prob, 
                   dt, 
                   method=method, 
                   adaptive=adaptive, 
                   tstops=tstops, 
                   always_new=always_new,
    )

    return sol
end

function _flat_reinit()
end

function _reinit(u0, p, sys)
    u_newton = copy(u0)
    u_homo = copy(u0)
    u_adap_homo = copy(u0)
    always_new = true
    solver_stats = false

    address = p[1]

    d_idx = address["balance_d"][sys.models.fault.bus[1]]
    q_idx = address["balance_q"][sys.models.fault.bus[1]]



    
    r1 = solve_newton!(u_newton, 
                       p, 
                       p[1]; 
                       always_new=always_new, 
                       solver_stats=solver_stats,
                       )
    r2 = solve_homotopy!(u_homo, 
                        p[1:end-1], 
                        p[1];
                        λ_target=1.0, 
                        Δλ=0.01, 
                        always_new=always_new,
                        vd_idx=d_idx,
                        vq_idx=q_idx,
                        )
    r3 = solve_adaptive_homotopy!(u_adap_homo, 
                                  p[1:end-1], 
                                  p[1]; 
                                  λ_target=1.0, 
                                  always_new=always_new,
                                  vd_idx=d_idx,
                                  vq_idx=q_idx,
                                )
    return (;r1, r2, r3, u_newton, u_adap_homo, u_homo)
end


function run_simulation(; case_bus=14, 
                load_factor=1.0, 
                method=:Euler,
                dt_list=[1e-5, 1e-6],
                verbose=false,
    )
    data_file = _get_data_file(case_bus)
    println("case file loaded.")
    sys = _create_system(data_file, load_factor)
    solve_power_flow!(sys)
    println("power flow solved.")
    run_static_init!(sys)
    _add_fault!(sys, case_bus)
    u0, p, mass_matrix = _create_dynamic_data(sys)
    sol_pf = _run_sim(u0, p, mass_matrix, t_end=0.001)
    println("steady state simulation done.")
    u_pf = sol_pf.u[end]
    println("entering dynamic simulation")
    # fault-on simulation
    sol_list=[]
    num_fails=0
    lambda = 1.0
    p = (p[1:end-1]..., lambda)
    for dt in dt_list
        sol = _run_sim(u_pf, p, mass_matrix, dt=dt, method=method, t_end=2*dt)
        if sol.retcode == :MaxIter
            num_fails += 1
        end
        push!(sol_list, sol)
    end
    println("Direct DAE integration failures = $num_fails")
    @info "Direct DAE integration failures = $num_fails"

    println("enterig re-init")
    
    ## flat re-init
    u_flat = copy(u_pf)
    u_flat[p[1]["balance_d"]] .= 1.0
    u_flat[p[1]["balance_q"]] .= 0.0
    # dt_flat = 1e-4
    dt_flat = 5e-4
    sol_flat = _run_sim(u_flat, p, mass_matrix, dt=dt_flat, method=method, t_end=dt_flat)
    
    # homotopy reinit
    reinit = _reinit(u_pf, p, sys)

    # post-re-init simulation
    lambda = 1.0
    p = (p[1:end-1]..., lambda)
    dt_post = 5e-4
    u0_post = copy(reinit.u_newton)
    println("entering post re-init run")
    # sol_post = _run_sim(u0_post, p, mass_matrix, dt=dt_post, method=method, t_end=dt_post)
    # sol_post = _run_sim(u0_post, p, mass_matrix, dt=dt_post, method=method, t_end=0.1)
    sol_post = _run_sim(u0_post, p, mass_matrix, dt=dt_post, method=method, t_end=0.0095)


    address = p[1]
    return (; sol_list, sol_post, reinit, sol_flat, sol_pf, address)
end
end #module