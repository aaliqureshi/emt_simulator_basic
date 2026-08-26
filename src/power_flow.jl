module PowerFlow

export solve_power_flow!

include("Models.jl"); using .Models

using NonlinearSolve, MyDiffEq

function power_flow_traditional!(du, u, p)
    T = eltype(u)

    address, models, Y, non_slack_buses, v_update_buses = p

    bus = models.bus
    generator = models.generator
    load = models.load

    V = Vector{T}(bus.v)
    theta = Vector{T}(bus.theta)

    V[v_update_buses] = @. u[address["balance_reactive_power"]]
    theta[non_slack_buses] = @. u[address["balance_real_power"]]

    P = zeros(T, length(bus.idx))
    Q = zeros(T, length(bus.idx))

    v = V .* cis.(theta)
    S = v .* conj.(Y * v)

    P = real.(S)
    Q = imag.(S)

    ## power balance : (power injection into network) + (load) - (generation) = 0
    P[generator.bus] -= @. generator.p_m
    P[load.bus] += @. load.p
    Q[load.bus] += @. load.q

    du[address["balance_reactive_power"]] .= @. Q[v_update_buses]
    du[address["balance_real_power"]] .= @. P[non_slack_buses]
end

function powerflow(u, p)
    T = eltype(u)
    du = zeros(T, length(u))
    power_flow_traditional!(du, u, p)
    return du
end

function solve_power_flow!(sys)
    # Assumes bus numbering is 1:n_bus (sys.Y and bus vectors indexed by bus number).
    models = sys.models
    Y = sys.Y
    non_slack_buses = sys.non_slack_buses
    v_update_buses = sys.v_update_buses

    # build address dict for power flow variables
    idx_balance_reactive_power = 1 : length(v_update_buses)
    idx_balance_real_power = idx_balance_reactive_power[end]+1 : idx_balance_reactive_power[end]+(length(non_slack_buses))

    address = Dict(
        "balance_reactive_power" => idx_balance_reactive_power,
        "balance_real_power" => idx_balance_real_power
    )

    # initial guess from bus data
    # u0 = vcat(models.bus.v[v_update_buses], models.bus.theta[non_slack_buses])
    #flat start
    u0 = vcat(ones(Float64, length(v_update_buses)), zeros(Float64, length(non_slack_buses)))
    # random start
    # a1=0.9
    # a2=1.1
    # b1=-2*pi*50/360
    # b2=2*pi*10/360
    # r1=a1.+(a2-a1).*rand(length(v_update_buses))
    # r2=b1.+(b2-b1).*rand(length(non_slack_buses))
    # u0=vcat(r1, r2)
    

    p = (address, models, Y, non_slack_buses, v_update_buses)

    # prob = MyDiffEq.NonLinearProblem(powerflow, u0, p)
    # sol = MyDiffEq.NLsolve(prob, method=:GradientDescent, 
    #                        atol = 1e-9, rtol = 1e-7, max_iter = 20000,alpha=1e-2)
    # u0 = sol.u_final


    # prob = NonlinearProblem(powerflow, u0, p)
    # sol = solve(prob, NewtonRaphson(); abstol = 1e-9, reltol = 1e-9, maxiters = 5000)
    
    # prob = MyDiffEq.NonLinearProblem(powerflow, u0, p)
    # sol = MyDiffEq.NLsolve(prob, method=:GradientDescent, 
    #                        atol = 1e-9, rtol = 1e-7, max_iter = 3500,alpha=1e-9)
    
    prob = MyDiffEq.NonLinearProblem(powerflow, u0, p)
    # sol = MyDiffEq.NLsolve(prob, 
    #                        method=:SRD,
    #                        max_iter = 75000,
    #                        atol = 1e-7,
    #                        rtol = 1e-7, 
    #                        alpha=1e-6,
    #                     #    step_decay=0.1,
    #                        verbose=true,
    #                     #    sampling_mode=:with_replacement,
    #                        sampling_mode=:random_reshuffling,
    #                     #    normalize=:kaczmarz,
    #                        )
    sol = MyDiffEq.NLsolve(prob, 
                        #    method=:LiftedGD,
                        # method=:GradientDescent,
                        method=:NewtonRaphson,
                           max_iter = 6000,
                        # max_iter=150,
                           atol = 1e-3,
                           rtol = 1e-5, 
                           alpha=1e-0,
                        # alpha = 1e-1,
                        # #    step_decay=0.1,
                           verbose=true,
                        # #    sampling_mode=:with_replacement,
                        #    sampling_mode=:random_reshuffling,
                        # #    normalize=:kaczmarz,
                        store_stats=false,
                           )
    
    # u0 = sol.u_final
    # @show sol.u_final
    # prob = NonlinearProblem(powerflow, u0, p)
    # sol = solve(prob, NewtonRaphson(); abstol = 1e-9, reltol = 1e-9, maxiters = 5000)
                       
    # populate models with the solution
    # models.bus.v[v_update_buses] .= sol.u[address["balance_reactive_power"]]
    # models.bus.theta[non_slack_buses] .= sol.u[address["balance_real_power"]]
    models.bus.v[v_update_buses] .= sol.u_final[address["balance_reactive_power"]]
    models.bus.theta[non_slack_buses] .= sol.u_final[address["balance_real_power"]]

    # convert to d-q components
    phasor2DP!(models.bus)

    return sol
end

end # module
