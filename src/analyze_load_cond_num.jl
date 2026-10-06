using Pkg; Pkg.activate(".")
using Revise
using Barq
using MyDiffEq, Plots

# 1. Load data
# data_file = "cases/Fault_Cases/ieee14_fault_barq_no_shunt.xlsx"
data_file = "cases/Fault_Cases/ieee39_fault.xlsx"
# data_file = "cases/Fault_Cases/SMIB_RL_Line_DrCui.xlsx"
# data_file = "cases/Simple_Cases/wecc_full_slack.xlsx"
# data_file = "cases/Fault_Cases/case2383wp_gc.xlsx"
# data_file = "cases/Fault_Cases/case118_gc.xlsx"
# data_file = "cases/Fault_Cases/case3012wp_barq.xlsx"
models = load_data(data_file)

# build system
sys = build_system(models)

# models.fault.bus = [20]
models.fault.bus = [24]
# models.fault.x_fault[1] = 0.006
models.fault.x_fault[1] = 0.006

# models.fault.bus = [1396]
# models.fault.x_fault[1] = 0.01

# case 118 fault
# models.fault.bus=[12]
# models.fault.x_fault[1]=0.009

# models.load.p[:] .*= 1.42

# models.load.p[end-553] = 4.0

# 2. Solve power flow
solve_power_flow!(sys);

# 3. Static initialization
run_static_init!(sys)

# 4. Dynamic simulation setup
address = build_dynamic_address(sys);
mass_matrix = build_mass_matrix(sys, address)
u0 = build_initial_conditions(sys, address);


function run_simulation(u0, lambda, dt; 
                        method=:Euler, 
                        always_new=true, 
                        tstops=[], 
                        adaptive=false, 
                        t_end=0.1,
    )
    p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)
    prob = MyDiffEq.ODEProblem(solve_dynamic_sim!, u0, (0.0, t_end), p, mass_matrix)
    sol = MyDiffEq.Solve(prob, dt, method=method, adaptive=adaptive, tstops=tstops, always_new=always_new)
    return sol
end

# Pre-fault simulation
lambda=0.0
dt0=5e-4
# method=:Trap
method=:Euler
sol_pf=run_simulation(u0, lambda, dt0, t_end=2*dt0, method=method);
u1 = sol_pf.u[end]

# fault-on simulation 3
dt_list = [5e-4, 5e-5, 5e-6]
# dt_list = [5e-4, 5e-5]
# dt_list = [5e-5, 5e-6]
# dt_list=[5e-6]
# dt_list = [1e-5, 1e-6]
# dt_list = [1e-4, 1e-5]
# dt_list = [1e-5, 5e-6]


lambda = 1.0
sol_list = []
num_fails = 0
# method=:Trap
for iter in eachindex(dt_list)
    ux = run_simulation(u1, lambda, dt_list[iter], t_end=dt_list[iter],method=method)
    # ux = run_simulation(u1, lambda, dt_list[iter], t_end=0.01,method=method)
    if ux.retcode == :MaxIter
        num_fails+=1
    end
    push!(sol_list, ux)
    # println("For dt=$(dt_list[iter]), retcode = $(ux.retcode)")
end
@show num_fails

using JSON

begin
    case="39"

    load_model = "P"
    # load_model = "I"
    # load_model = "Z"
    # load_model = "ZIP"

    # PQBRAK="N"
    PQBRAK="Y"

    con = "y"

    filename="$(case)_$(load_model)_$(PQBRAK).json"

    # filename="$(case)_$(load_model)_$(PQBRAK)_$(con)"
end

open(filename, "w") do io
    JSON.print(io, sol_list[2].newton_log)
end







p = open("39_P_N.json", "r") do io
    JSON.parse(io)
end
z = open("39_Z_N.json", "r") do io
    JSON.parse(io)
end
i = open("39_I_N.json", "r") do io
    JSON.parse(io)
end
zip = open("39_ZIP_N.json", "r") do io
    JSON.parse(io)
end


p_brak = open("39_P_Y.json", "r") do io
    JSON.parse(io)
end
z_brak = open("39_Z_Y.json", "r") do io
    JSON.parse(io)
end
i_brak = open("39_I_Y.json", "r") do io
    JSON.parse(io)
end
zip_brak = open("39_ZIP_Y.json", "r") do io
    JSON.parse(io)
end
using Plots

plot(p.cond[1:10], yscale=:log10)
plot!(z.cond, yscale=:log10)
plot!(i.cond, yscale=:log10)
plot!(zip.cond[1:10], yscale=:log10)

plot(p_brak.cond[1:10], yscale=:log10)
plot!(z_brak.cond, yscale=:log10)
plot!(i_brak.cond, yscale=:log10)
plot!(zip_brak.cond[1:10], yscale=:log10)

plot(p.cond[1:10], yscale=:log10, label="NO PQBRAK", linewidth=2.5, title="P")
plot!(p_brak.cond[1:10], yscale=:log10, label="PQBRAK", linewidth=2.5,)

plot(z.cond, yscale=:log10, label="NO PQBRAK", linewidth=2.5, title="Z")
plot!(z_brak.cond, yscale=:log10, label="PQBRAK", linewidth=2.5,)

plot(i.cond, yscale=:log10, label="NO PQBRAK", linewidth=2.5, title="I")
plot!(i_brak.cond, yscale=:log10, label="PQBRAK", linewidth=2.5,)

plot(zip.cond, yscale=:log10, label="NO PQBRAK", linewidth=2.5, title="ZIP")
plot!(zip_brak.cond, yscale=:log10, label="PQBRAK", linewidth=2.5,)




## algeb cond
plot(p.cond_algeb, yscale=:log10, label="NO PQBRAK", linewidth=2.5, title="P")
plot!(p_brak.cond_algeb, yscale=:log10, label="PQBRAK", linewidth=2.5,)

plot(z.cond_algeb, yscale=:log10, label="NO PQBRAK", linewidth=2.5, title="Z")
plot!(z_brak.cond_algeb, yscale=:log10, label="PQBRAK", linewidth=2.5,)

plot(i.cond_algeb, yscale=:log10, label="NO PQBRAK", linewidth=2.5, title="I")
plot!(i_brak.cond_algeb, yscale=:log10, label="PQBRAK", linewidth=2.5,)

plot(zip.cond_algeb, yscale=:log10, label="NO PQBRAK", linewidth=2.5, title="ZIP")
plot!(zip_brak.cond_algeb, yscale=:log10, label="PQBRAK", linewidth=2.5,)





plot(zip.cond, yscale=:log10, label="DAE cond", linewidth=2.5, title="ZIP - NO PQBRAK")
plot!(zip.cond_algeb, yscale=:log10, label="Algeb cond", linewidth=2.5, title="ZIP - NO PQBRAK")

plot(zip_brak.cond, yscale=:log10, label="DAE cond", linewidth=2.5, title="ZIP - PQBRAK")
plot!(zip_brak.cond_algeb, yscale=:log10, label="Algeb cond", linewidth=2.5, title="ZIP - PQBRAK")


plot(zip.sigma_min, yscale=:log10, label="NO PQBRAK", linewidth=2.5, title="ZIP - sigma min")
plot!(zip_brak.sigma_min, yscale=:log10, label="PQBRAK", linewidth=2.5,)

plot(zip.sigma_max, yscale=:log10, label="NO PQBRAK", linewidth=2.5, title="ZIP - sigma max")
plot!(zip_brak.sigma_max, yscale=:log10, label="PQBRAK", linewidth=2.5,)



lambda = 1.0
sol_list = []
num_fails = 0
# method=:Trap
sol_brak = run_simulation(u1, lambda, 5e-4, t_end=5e-4,method=method)
sol_zip = run_simulation(u1, lambda, 5e-4, t_end=5e-4,method=method)

using LinearAlgebra

brak_zip_diff = norm(sol_brak.u[end] - sol_zip.u[end], Inf)

brak_zip_diff = norm(sol_brak.u[end] - sol_zip.u[end], 2)

@show sol_zip.u[end][address["balance_d"]]

@show sol_brak.u[end][address["balance_d"]]

norm(sol_zip.u[end][address["balance_d"]] - sol_brak.u[end][address["balance_d"]], Inf)

norm(sol_zip.u[end][address["fault_id"]] - sol_brak.u[end][address["fault_id"]], Inf)

norm(sol_zip.u[end][address["fault_iq"]] - sol_brak.u[end][address["fault_iq"]], Inf)


sol = run_simulation(u1, lambda, 5e-4, t_end=5e-4,method=method)



zip_brak1=copy(sol.newton_log.cond_algeb)

plot(zip.cond_algeb, yscale=:log10, label="NO PQBRAK", linewidth=2.5, title="ZIP")
plot!(zip_brak.cond_algeb, yscale=:log10, label="PQBRAK 0.7", linewidth=2.5,)
plot!(zip_brak6, yscale=:log10, label="PQBRAK 0.6", linewidth=2.5,)
plot!(zip_brak5, yscale=:log10, label="PQBRAK 0.5", linewidth=2.5,)
plot!(zip_brak1, yscale=:log10, label="PQBRAK 0.1", linewidth=2.5,)
