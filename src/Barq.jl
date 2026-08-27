module Barq

include("Models.jl");      using .Models
include("utils.jl");       using .Utils
include("data_loader.jl"); using .DataLoader
include("system.jl");      using .SystemModel
include("power_flow.jl");  using .PowerFlow
include("static_init.jl"); using .StaticInit
include("dynamic_sim.jl"); using .DynamicSim
include("homotopy.jl");    using .Homotopy

# device models
export Bus, Line, Generator, Fault, Load, Slack
export solve_generator!, solve_line!, solve_fault!, balance!
export phasor2DP!, compute_line_currents!, compute_load_currents!

# network construction
export load_data, System, build_system
export build_Y_matrix, build_incidence_matrix, build_B_matrix

# initialization and simulation
export solve_power_flow!, run_static_init!
export solve_dynamic_sim!, build_dynamic_address, build_mass_matrix, build_initial_conditions

# algebraic solvers and continuation
export solve_newton!, solve_damped_newton!, solve_backtracking_newton!
export solve_levenberg_marquardt!, solve_homotopy!, solve_homotopy_lm!, solve_adaptive_homotopy!
export solve_algebraic!, pseudo_arclength!

end # module
