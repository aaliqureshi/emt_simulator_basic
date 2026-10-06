using ForwardDiff, LinearAlgebra, Printf
using OrdinaryDiffEq


# function _load_current(i_L0, v_d, v_q; v_min=0.7)
# end

# steady state GFM operation
function GFM!(du, u, p, t)
    e_q, i_ld, i_lq, v_d, v_q, i_cd, i_cq = u
    v_ref, v_slack, r_line, x_line, ω, tq, ki, kp, i_L, g_f, i_max = p
    
    V = sqrt(v_d^2 + v_q^2)
    
    # constnat current model with conversion to constant impedance when V<0.7
    i_Ld = i_L*v_d/(max(V, 0.7))
    i_Lq = i_L*v_q/(max(V, 0.7))

    # fault current
    i_fd = g_f*v_d
    i_fq = g_f*v_q

    du[1] = (ki*(v_ref - V) - e_q)/tq
    du[2] = (v_d - v_slack - r_line*i_ld + x_line*i_lq)*ω/x_line
    du[3] = (v_q - r_line*i_lq - x_line*i_ld)*ω/x_line
    du[4] = i_cd - i_ld - i_Ld - i_fd
    du[5] = i_cq - i_lq - i_Lq - i_fq
    du[6] = v_d - v_ref
    du[7] = v_q
    
    # return du
end

# current-limited model
function GFL!(du, u, p, t)
    e_q, i_ld, i_lq, v_d, v_q, i_cd, i_cq = u
    v_ref, v_slack, r_line, x_line, ω, tq, ki, kp, i_L, g_f, i_max = p

    V = sqrt(v_d^2 + v_q^2)

    # constnat current model with conversion to constant impedance when V<0.7
    i_Ld = i_L*v_d/(max(V, 0.7))
    i_Lq = i_L*v_q/(max(V, 0.7))

    # fault current
    i_fd = g_f*v_d
    i_fq = g_f*v_q

    du[1] = (ki*(v_ref - V) - e_q)/tq
    du[2] = (v_d - v_slack - r_line*i_ld + x_line*i_lq)*ω/x_line
    du[3] = (v_q - r_line*i_lq - x_line*i_ld)*ω/x_line
    du[4] = i_cd - i_ld - i_Ld - i_fd
    du[5] = i_cq - i_lq - i_Lq - i_fq
    du[6] = i_cq - (e_q + kp*(v_ref - V))
    du[7] = i_cd^2 + i_cq^2 - i_max^2
    
end


function Algeb_Reinit!(du, u, p, t)
    i_ld, i_lq, v_d, v_q, i_cd, i_cq = u
    v_ref, v_slack, r_line, x_line, ω, tq, ki, kp, i_L, g_f, i_max = p
 
    e_q = 0.0 # steady state of e_q

    V = sqrt(v_d^2 + v_q^2)
    
    # constnat current model with conversion to constant impedance when V<0.7
    i_Ld = i_L*v_d/(max(V, 0.7))
    i_Lq = i_L*v_q/(max(V, 0.7))

    # fault current
    i_fd = g_f*v_d
    i_fq = g_f*v_q

    du[1] = (v_d - v_slack - r_line*i_ld + x_line*i_lq)*ω/x_line
    du[2] = (v_q - r_line*i_lq - x_line*i_ld)*ω/x_line
    du[3] = i_cd - i_ld - i_Ld - i_fd
    du[4] = i_cq - i_lq - i_Lq - i_fq
    du[5] = i_cq - (e_q + kp*(v_ref - V))
    du[6] = i_cd^2 + i_cq^2 - i_max^2
    
end


# parameters v_ref, v_slack, r_line, x_line, ω, tq, ki, kp, i_L, g_f, i_max
# p = (1.0, 1.0, 0.08, 0.35, 2π*60, 0.2, 6.5, 4.7, 0.32, 0.6, 1.43)
p = (1.0, 1.0, 0.08, 0.35, 2π*60, 0.2, 6.5, 4.7, 0.32, 1.119, 1.43)
# p = (1.0, 1.0, 0.08, 0.35, 2π*60, 0.2, 6.5, 4.7, 0.32, 0.6, 0.9)
# p= (1.0, 1.0, 0.08, 0.35, 2π*60, 0.2, 6.5, 8.0, 0.2, 1.70, 1.60)
# p= (1.0, 1.0, 0.08, 0.35, 2π*60, 0.2, 6.5, 8.0, 0.2, 0.0, 1.60)

# g_f = 1.4889
# p = (1.0, 1.0, 0.08, 0.35, 2π*60, 0.2, 6.5, 8.0, 0.32, g_f, 1.43)
p_ss = (p[1:end-2]..., 0.0, p[end])
u0 = ones(7)
mass_matrix = zeros(7,7)
prob0_ss = OrdinaryDiffEq.ODEFunction(GFM!, mass_matrix=mass_matrix) 
prob1_ss = OrdinaryDiffEq.ODEProblem(prob0_ss, u0, (0.0, 0.1), p_ss)
sol0_ss = OrdinaryDiffEq.solve(prob1_ss, 
                      Trapezoid(), 
                      adaptive=false,
                      dt=1e-4,
                    #   initializealg = SciMLBase.NoInit(),
                      )

u1 = sol0_ss.u[end]

mass_matrix = zeros(7,7)
mass_matrix[1,1] = 1.0
# mass_matrix[2,2] = 1.0
# mass_matrix[3,3] = 1.0
prob0_d = OrdinaryDiffEq.ODEFunction(GFL!, mass_matrix=mass_matrix) 
prob1_d = OrdinaryDiffEq.ODEProblem(prob0_d, u1, (0.0, 0.1), p)
sol0_d = OrdinaryDiffEq.solve(prob1_d, 
                      Trapezoid(), 
                      adaptive=false,
                      dt=1e-4,
                      initializealg = SciMLBase.NoInit(),
                      )
# prob1 = MyDiffEq.ODEProblem(GFL!, u1, (0.0, 0.01), p, mass_matrix)
# prob1 = MyDiffEq.ODEProblem(GFM!, u1, (0.0, 0.01), p, mass_matrix)
# sol1 = MyDiffEq.Solve(prob1, 
#                     0.0001, 
#                     method=:Euler, 
#                     adaptive=false,
#                     tstops=[],
#                     always_new=true,
#                     )
@show sol1.u[end]

# dt_list = [1e-4, 1e-5, 1e-6]
# dt_list = [5e-4, 5e-5, 5e-6]
# dt_list = [1e-5, 1e-6]
# dt_list = [5e-3, 5e-4, 5e-5, 5e-6]
# dt_list = [5e-3, 5e-4, 5e-5, 1e-6]
sol_list = []
num_fails = 0
method = :Euler
# prob2 = MyDiffEq.ODEProblem(GFL!, u1, (0.0, 0.5), p, mass_matrix)
for iter in eachindex(dt_list)
    # prob2 = MyDiffEq.ODEProblem(GFL!, u1, (0.0, dt_list[iter]), p, mass_matrix)
    prob0_d = OrdinaryDiffEq.ODEFunction(GFL!, mass_matrix=mass_matrix) 
    prob1_d = OrdinaryDiffEq.ODEProblem(prob0_d, u1, (0.0, 0.1), p)
    ux = OrdinaryDiffEq.solve(prob1_d, 
                          Trapezoid(), 
                          adaptive=false,
                          dt=dt_list[iter],
                          initializealg = SciMLBase.NoInit()
                          )
    if (ux.retcode) != 1
        println("failed at dt = $(dt_list[iter])")
    end
    
    push!(sol_list, ux)
    # println("For dt=$(dt_list[iter]), retcode = $(ux.retcode)")
end
num_fails

nr_res_norm = []
#local implementation of NR method
function solve_nr(residual_func!, initial_guess, p; max_iter=20, tol=1e-8, verbose=true)
    u = float.(collect(initial_guess))
    r = zeros(length(u))
    for k in 1:max_iter
        residual_func!(r, u, p, 0)
        push!(nr_res_norm, norm(r, 2))
        any(!isfinite, r) && return nothing
        J = ForwardDiff.jacobian(x -> begin
                                        du = similar(x)
                                        residual_func!(du, x, p, 0)
                                        du
                                      end, u)
        verbose && @printf("  iter %2d | residual norm: %.4e\n  | cond: %.4e", k, norm(r, 2), cond(J))
        norm(r, Inf) <= tol && return u, norm(r, 2)

        # J = ForwardDiff.jacobian(x -> begin
        #                                 du = similar(x)
        #                                 residual_func!(du, x, p, 0)
        #                                 du
        #                               end, u)
        d = J \ (-r)
        any(!isfinite, d) && return norm(r, 2)
        u += d
    end
    residual_func!(r, u, p, 0)
    return norm(r, Inf) <= tol ? u : nr_res_norm
end

"""Parameter tuple with the fault conductance g_f replaced by gf."""
_with_gf(p, gf) = (p[1:9]..., gf, p[11])

res_hist = []
lambda_hist = []
u_hist = []
"""
Natural-parameter continuation from the unfaulted to the faulted network: the fault
conductance is ramped as g_f(lambda) = lambda*g_f, so lambda enters the Jacobian.
Returns nothing if any continuation step fails.
"""
function solve_homotopy(residual_func!, u0, p; nsteps=5, kwargs...)
    u = float.(collect(u0))
    # for lambda in range(0.0, 1.0, length=nsteps+1)
    icd = Float64[]
    icq = Float64[]
    for lambda in collect(0.0:0.01:1.0)
        push!(lambda_hist, lambda)
        unew, res_norm = solve_nr(residual_func!, u, _with_gf(p, lambda*p[10]); kwargs...)
        # @show unew
        if isnothing(unew)
            @printf("Homotopy failed at lambda = %.3f\n", lambda)
            return nothing
        end
        @printf("lambda = %.3f  ->  V = %.6f, |i_c| = %.6f\n",
                lambda, hypot(unew[3], unew[4]), hypot(unew[5], unew[6]))
        u = unew
        push!(u_hist, u)
        push!(res_hist, res_norm)
        push!(icd, u[5])
        push!(icq, u[6])
    end
    return u, icd, icq
end

u2 = u1[2:end]
nr_res_norm = solve_nr(Algeb_Reinit!, u2, p)

kk_homo, icd, icq = solve_homotopy(Algeb_Reinit!, u2, p; verbose=false);



# dt_list = [1e-4, 1e-5, 1e-6]
dt = 1e-4
prob0_homo = OrdinaryDiffEq.ODEFunction(
                GFL!, 
                mass_matrix=mass_matrix,
                )
prob1_homo = OrdinaryDiffEq.ODEProblem(prob0_homo,
                            vcat(u2[1], kk_homo), 
                            (0.0, 0.5), 
                            p,
                            )
sol3 = OrdinaryDiffEq.solve(prob1_homo,
                    Trapezoid(),
                    dt=dt, 
                    # method=:Euler, 
                    adaptive=false,
                    # tstops=[],
                    # always_new=true,
                    initializealg = SciMLBase.NoInit(),
                    )


using Plots, LaTeXStrings

default(
    fontfamily = "Computer Modern",
    linewidth = 3.5,
    markersize = 7,
    markerstrokewidth = 0,
    legendfontsize = 13,
    guidefontsize = 14,
    tickfontsize = 14,
    titlefontsize = 14,
    dpi = 300,
    # size = (1000, 380)
)
e=15
e2 = length(sol3.newton_log.residual_norm)
pa = plot(
    1:e,
    sol_list[1].newton_log.residual_norm[1:e],
    yscale = :log10,
    color = :blue,
    linestyle = :dash,
    marker = :diamond,
    # label = "No re-init. (h = 500 μs)",
    label = "No re-init. (h = $(trunc(Int, dt_list[1]/1e-6)) μs)"
)
plot!(
    pa,
    1:e,
    sol_list[2].newton_log.residual_norm[1:e],
    yscale = :log10,
    color = :red,
    linestyle = :dash,
    marker = :diamond,
    # label = "No re-init. (h = 500 μs)",
    label = "No re-init. (h = $(trunc(Int, dt_list[2]/1e-6)) μs)"
)
plot!(
    pa,
    1:e,
    sol_list[3].newton_log.residual_norm[1:e],
    yscale = :log10,
    color = :green,
    linestyle = :dash,
    marker = :diamond,
    # label = "No re-init. (h = 500 μs)",
    label = "No re-init. (h = $(trunc(Int, dt_list[3]/1e-6)) μs)"
)

plot!(
    pa,
    1:e2,
    sol3.newton_log.residual_norm[1:e2],
    color = :teal,
    linestyle = :solid,
    marker = :square,
    label = "With re-init.",
)

plot!(
    pa,
    xlabel = "Newton iteration",
    ylabel = L"\|F(\mathbf{z})\|_2",
    framestyle = :box,
    grid = true,
    gridalpha = 0.18,
    minorgrid = false,
    legend = :bottomright,
    # legend = (0.42,0.42),
    title = "Post-event Integration (implicit $method)",
    # titlelocation = :left,
    # xlims = (0.7, xmax_a + 0.3),   # small padding on both sides
    # left_margin = 6Plots.mm,
    # bottom_margin = 5Plots.mm,
)



# NR Re-init norm
e3 = 15
pb = plot(
    1:e3,
    nr_res_norm[1:e],
    yscale = :log10,
    color = :blue,
    linestyle = :dash,
    marker = :diamond,
    # label = "No re-init. (h = 500 μs)",
    # label = "NR-based re-init."
    label="",
)
plot!(
    pb,
    xlabel = "Newton iteration",
    ylabel = L"\|G(\mathbf{z})\|_2",
    framestyle = :box,
    grid = true,
    gridalpha = 0.18,
    minorgrid = false,
    legend = :bottomright,
    # legend = (0.42,0.42),
    title = "NR-based re-init.",
    # titlelocation = :left,
    # xlims = (0.7, xmax_a + 0.3),   # small padding on both sides
    # left_margin = 6Plots.mm,
    # bottom_margin = 5Plots.mm,
)