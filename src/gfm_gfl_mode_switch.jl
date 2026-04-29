using Printf, Plots

# --- System Parameters ---
const V_ref  = 1.0  # GFM internal voltage reference (pu)
const R_g    = 0.1  # Output virtual/physical impedance (pu)
const P_load = 0.9  # Constant power demand (pu)
const I_max  = 3.0  # Saturated fault current limit (GFL mode)

# --- Manifold 1: Grid Forming Mode ---
function gfm_residual(V)
    return (V_ref - V) / R_g - (P_load / V)
end

function gfm_jacobian(V)
    return -1.0 / R_g + (P_load / (V^2))
end

function gfl_residual(V)
   return I_max - (P_load / V)
end

function gfl_jacobian(V)
    return (P_load / V^2)
end

function homotopy_residual_mode_switch(V, lambda)
    return (1-lambda) * gfm_residual(V) + lambda * gfl_residual(V)
end

function homotopy_jacobian_mode_switch(V, lambda)
    return (1-lambda)*gfm_jacobian(V) + lambda*gfl_jacobian(V)
end

# --- Newton-Raphson Solver ---
function solve_nr(residual_func, jacobian_func, initial_guess; max_iter=10)
    V = initial_guess
    @printf("Initial Guess: V = %.4f pu\n", V)
    
    for i in 1:max_iter
        res = residual_func(V)
        jac = jacobian_func(V)
        
        @printf("Iter %d: V = %8.4f | Residual = %8.4f | Jacobian = %8.4f\n", i, V, res, jac)
        
        if abs(res) < 1e-6
            println("--> Converged!\n")
            return V
        end
        
        # The Convergence Killer
        if abs(jac) < 1e-12
            println("--> FATAL ERROR: Jacobian is singular (divide by zero)!\n")
            return NaN
        end
        
        # Newton Update Step
        V = V - (res / jac)
    end
    println("--> Failed to converge after max iterations.\n")
    return V
end

# --- Robust NR Solver for Homotopy Steps ---
function solve_nr_homotopy(V_guess, lambda; max_iter=10)
    V = V_guess
    for i in 1:max_iter
        res = homotopy_residual_mode_switch(V, lambda)
        jac = homotopy_jacobian_mode_switch(V, lambda)
        
        if abs(res) < 1e-8
            println("--> $(lambda) ---> Converged!\n")
            return V
        end
        V = V - (res / jac)
    end
    return NaN # Failure
end




####
lw = 2.5
ms = 7.0
V = 0.1:0.01:1.0
# V = -0.0:0.01:3.0
vv = @. gfm_residual(V)
gfm_curve= vv
kk = @. 3.0 - P_load/V
pa = plot(V, vv, lw=lw, label= "GFM Mode")
plot!(V, kk, lw=lw, label= "GFL mode - Current residual")
scatter!([0.9],[0.0], marker=:diamond, markersize=ms, markerstrokewidth=0, label="pre-event solution")
scatter!([0.3], [0.0], marker=:diamond, markersize=ms, markerstrokewidth=0, label="post-event solution")
scatter!([0.9], [gfl_residual(0.9)], marker=:diamond, markersize=ms, markerstrokewidth=0, label="pre-event guess on post-event manifold")
plot!(pa, xlabel="Voltage (p.u)",
      ylabel="g(z)",
      legendfontsize=6,
         tickfontsize=8, guidefontsize=9,
         framestyle=:box, grid=true, gridalpha=0.2)
# fig = plot(pa, layout=(1, 1), size=(700, 300), dpi=300,
#          margin=4Plots.mm, bottom_margin=5Plots.mm)
fig = plot(pa, dpi=300)
savefig(fig, "converter_mode_switch.pdf")
####

###
# lambdas = 0.0:0.1:1.0
# cc = @. homotopy_residual_mode_switch(0.9, lambdas)
###

# ==========================================
# SIMULATION SCRIPT
# ==========================================

println("=== Phase 1: Pre-Fault (GFM Mode) ===")
# Starting simulation from cold start with a good guess (V=1.0)
V_prefault = solve_nr(gfm_residual, gfm_jacobian, 1.0)

println("=== Phase 2: During Fault (GFL Mode) ===")
# Network topology changes or converter hits limits. 
# It switches to the GFL manifold (Constant Current).
# We analytically find the state on this new manifold: I_max - P/V = 0

# V_fault = P_load / I_max
V_fault = solve_nr(gfl_residual, gfl_jacobian, V_prefault)
# @printf("System settles on GFL manifold at V = %.4f pu\n\n", V_fault)

# Number of continuation steps
n_steps = 10
# lambdas = range(0.0, stop=1.0, length=n_steps+1)
lambdas = 0.0:0.1:1.0

V0 = V_prefault
success = true
for (i, λ) in enumerate(lambdas)
    global V0
    V_next = solve_nr_homotopy(V0, λ)
    
    if isnan(V_next)
        @printf("Step %d (lambda=%.2f): FAILED to converge!\n", i-1, λ)
        success = false
        break
    else
        @printf("Step %d (lambda=%.2f): Converged to V = %.4f\n", i-1, λ, V_next)
        V0 = V_next # Update the guess for the NEXT lambda step
    end
end

println("=== Phase 3: Fault Cleared (Switching back to GFM Manifold) ===")
# The fault clears. We swap the algebraic equations back to GFM.
# The solver has no choice but to use the last known state (V_fault) as its initial guess.
V_postfault = solve_nr(gfm_residual, gfm_jacobian, V_fault)



######
using Printf

# --- System Parameters (Same as before) ---
const V_ref  = 1.0
const R_g    = 0.1
const P_load = 0.9
const I_max  = 3.0

# --- Manifold Definitions ---
g_fault(V)  = I_max - (P_load / V)
g_target(V) = (V_ref - V) / R_g - (P_load / V)

# --- Homotopy Function and its Jacobian ---
# H(V, lambda) = (1-lambda)*g_fault + lambda*g_target
function homotopy_residual(V, lambda)
    return (1 - lambda) * g_fault(V) + lambda * g_target(V)
end

function homotopy_jacobian(V, lambda)
    # Derivative of H with respect to V
    j_fault  = P_load / (V^2)
    j_target = -1.0 / R_g + (P_load / (V^2))
    return (1 - lambda) * j_fault + lambda * j_target
end

# --- Robust NR Solver for Homotopy Steps ---
function solve_nr_step(V_guess, lambda; max_iter=10)
    V = V_guess
    for i in 1:max_iter
        res = homotopy_residual(V, lambda)
        jac = homotopy_jacobian(V, lambda)
        
        if abs(res) < 1e-8
            return V
        end
        V = V - (res / jac)
    end
    return NaN # Failure
end

# ==========================================
# HOMOTOPY EXECUTION
# ==========================================

println("=== Transitioning from GFL to GFM via Homotopy ===")
# Starting point: The saturated fault voltage
V_current = P_load / I_max 
@printf("Starting Point (lambda=0.0): V = %.4f\n\n", V_current)

# Number of continuation steps
n_steps = 10
lambdas = range(0.0, stop=1.0, length=n_steps+1)

success = true
for (i, λ) in enumerate(lambdas)
    global V_current
    V_next = solve_nr_step(V_current, λ)
    
    if isnan(V_next)
        @printf("Step %d (lambda=%.2f): FAILED to converge!\n", i-1, λ)
        success = false
        break
    else
        @printf("Step %d (lambda=%.2f): Converged to V = %.4f\n", i-1, λ, V_next)
        V_current = V_next # Update the guess for the NEXT lambda step
    end
end

if success
    println("\n[SUCCESS] Homotopy successfully bridged the manifold gap!")
    @printf("Final Post-Fault Equilibrium: V = %.4f pu\n", V_current)
else
    println("\n[FAILURE] Even homotopy couldn't save this step size.")
end




####
using Printf
using Plots
using LaTeXStrings

# Set global plotting style for IEEE publication standard
default(
    fontfamily = "Computer Modern", # LaTeX font
    # titlesize = 12,
    guidefont = 10,
    tickfont = 8,
    legendfont = 8,
    grid = true,
    framestyle = :box,
    thickness_scaling = 1.1, # Make lines a bit thicker for clarity
    palette = :tab10 # Good distinct colors
)

# --- Define the System (Identical to MWE) ---
const V_ref  = 1.0  # GFM internal voltage reference (pu)
const R_g    = 0.1  # Output virtual/physical impedance (pu)
const P_load = 0.9  # Constant power demand (pu)
const I_max  = 3.0  # Saturated fault current limit (GFL mode)

# Pre-calculate crucial points
V_singular = 0.3
V_fault    = P_load / I_max # Also 0.3
# Analytical solutions of GFM: V^2 - V_ref*V + P*R = 0
V_stable   = 0.9 
V_unstable = 0.1

# define the range for voltage plotting
V_range = 0.05:0.005:1.1

# ------------------------------------------------------------------
# Visualization Functions
# ------------------------------------------------------------------

# Physical manifolds: solutions exist where I_inj = I_load
I_grid_inj(V) = (V_ref - V) / R_g # The GFM line
I_load_dem(V) = P_load / V        # The load curve

# Homotopy setup (needed to trace the path)
g_fault(V)  = I_max - (P_load / V)
g_target(V) = (V_ref - V) / R_g - (P_load / V)

function homotopy_residual(V, λ)
    return (1 - λ) * g_fault(V) + λ * g_target(V)
end

function homotopy_jacobian(V, λ)
    j_fault  = P_load / (V^2)
    j_target = -1.0 / R_g + (P_load / (V^2))
    return (1 - λ) * j_fault + λ * j_target
end

# Enhanced NR step solver that returns success and path data
function solve_nr_step_path(V_guess, λ)
    V = V_guess
    path_data = [V] # Store path during *this* step
    for i in 1:20 # Allow more steps inside homotopy
        res = homotopy_residual(V, λ)
        jac = homotopy_jacobian(V, λ)
        
        if abs(res) < 1e-8
            return (V, path_data, true)
        end
        # Safety for inner step (though homotopy should prevent this)
        if abs(jac) < 1e-12; return (NaN, path_data, false); end 
        
        V = V - (res / jac)
        push!(path_data, V)
    end
    return (NaN, path_data, false)
end


# ==================================================================
# Execute Simulation and Generate Plots
# ==================================================================

println("\n=== Generating IEEE-Style Figures ===\n")

# -- Initialize the Plot Layout (single column width) --
l = @layout [a ; b]
final_plot = plot(layout = l, size=(450, 700)) # 450 pixels ~ 3.25 inches

# ------------------------------------------------------------------
# PANEL A: The Newton-Raphson Singularity (V_fault = 0.3)
# ------------------------------------------------------------------
p1 = plot(title=L"(a) ~ Direct~NR~at~singular~post\text{-}fault~initial~guess", 
          xlabel=L"Converter~Terminal~Voltage,~V ~ \text{[pu]}", 
          ylabel=L"Current~ \text{[pu]}")

# 1. Plot background physical manifolds
plot!(p1, V_range, I_load_dem.(V_range), label=L"Load~Demand~(P/V)", color=:gray, style=:dash, alpha=0.5)
plot!(p1, V_range, I_grid_inj.(V_range), label=L"Grid~Injection~((V_{ref}-V)/R_g)", color=:black, style=:solid)

# 2. Mark the relevant equilibrium points
# Stable GFM solution (the target)
scatter!(p1, [V_stable], [I_load_dem(V_stable)], color=:forestgreen, marker=:circle, markersize=6, label=L"Target~Stable~Point")
# Initial condition (the GFL fault point)
scatter!(p1, [V_fault], [I_load_dem(V_fault)], color=:firebrick, marker=:xcross, markersize=8, label=L"Initial~Guess~(V_{fault})")

# 3. CRUCIAL: Show the tangent line (Jacobian) at the failure point
# The function NR solves is R(V) = I_inj - I_load. Solutions are where R(V)=0.
# The Jacobian is J = -1/R - (-P/V^2) = -10 + 0.9/0.3^2 = -10 + 0.9/0.09 = -10 + 10 = 0.
V_near_fault = V_fault-0.1 : 0.01 : V_fault+0.1
tangent_at_fault = fill(I_grid_inj(V_fault), length(V_near_fault))
plot!(p1, V_near_fault, tangent_at_fault, label=L"NR~Tangent~(J=0.0)", color=:gold, style=:dot, width=2.5)

# Annotation for clarity
annotate!(p1, V_fault+0.15, I_grid_inj(V_fault)-0.2, text(L"J(V_{fault})=0 \rightarrow Fail!", 8, :gold, :left))

# Formatting p1
xlims!(p1, (0.0, 1.1))
ylims!(p1, (0.0, 10.0))
# Move legend inside to save space
plot!(p1, legend=:topright)

# Add p1 to layout
plot!(final_plot[1], p1)
println("-> Panel A created.")


# ------------------------------------------------------------------
# PANEL B: Homotopy Path to Healthy Solution
# ------------------------------------------------------------------
p2 = plot(title=L"(b) ~ Homotopy~continuation~bridging~equilibrium~manifolds", 
          xlabel=L"Converter~Terminal~Voltage,~V ~ \text{[pu]}", 
          ylabel=L"Current~ \text{[pu]}")

# 1. Plot same background physical manifolds
plot!(p2, V_range, I_load_dem.(V_range), label="", color=:gray, style=:dash, alpha=0.5)
plot!(p2, V_range, I_grid_inj.(V_range), label="", color=:black, style=:solid)

# 2. Execute and Trace the Homotopy
n_steps = 10
lambdas = range(0.0, stop=1.0, length=n_steps+1)
V_current = V_fault # Start at 0.3

# Initialize storage for the *successful* continuous homotopy path
total_V_path = Float64[]
total_I_load_path = Float64[]
push!(total_V_path, V_current)
push!(total_I_load_path, I_load_dem(V_current))

success = true
println("-> Executing Homotopy...")
for (i, λ) in enumerate(lambdas)
    global V_current
    V_next, inner_path_V, inner_success = solve_nr_step_path(V_current, λ)
    
    if !inner_success
        println("   --> Step failed at lambda=", λ)
        success = false; break
    else
        @printf("   --> Step %d (lambda=%.2f): Converged to V = %.4f\n", i-1, λ, V_next)
        
        # Capture the whole inner solver path to make the line continuous
        append!(total_V_path, inner_path_V)
        for v in inner_path_V
            push!(total_I_load_path, I_load_dem(v))
        end
        
        V_current = V_next # Update start point for next lambda step
    end
end

# 3. Plot the trace of the Homotopy solution
# if success
    # Plot the full continuous path
    plot!(p2, total_V_path, total_I_load_path, label=L"Homotopy~Trace~(V \text{ vs } \lambda)", color=:royalblue, arrow=:arrow, width=2)
    
    # Trace arrows are good, but IEEE reviewers like marking the *successful* NR steps explicitly
    # Mark where λ=0.0, 0.1, 0.2 ... 1.0 landed on the curve. We can't use the simple linear lambdas 
    # to find the specific voltage point without re-running, so we just sample the continuous path.
    sampled_indices = round.(Int, range(1, stop=length(total_V_path), length=n_steps+1))
    scatter!(p2, total_V_path[sampled_indices], total_I_load_path[sampled_indices], color=:dodgerblue, marker=:diamond, markersize=3, label=L"NR~Sub\text{-}steps")

    # Re-mark start/end points explicitly
    scatter!(p2, [V_fault], [I_load_dem(V_fault)], color=:firebrick, marker=:xcross, markersize=8, label=L"Fault~State~(\lambda=0)")
    scatter!(p2, [V_stable], [I_load_dem(V_stable)], color=:forestgreen, marker=:circle, markersize=6, label=L"Target~Solution~(\lambda=1)")
    
    # Formatting p2
    # xlims!(p2, (0.0, 1.1))
    # ylims!(p2, (0.0, 10.0))
    plot!(p2, legend=:topright)
    
    # Add p2 to layout
    plot!(final_plot[2], p2)
    println("-> Panel B created.")
# else
    # println("[ERROR] Homotopy failed.")
# end


# ==================================================================
# Final Save for IEEE publication
# ==================================================================

# Save as PDF (or EPS) for high resolution vector graphics
# (requires pdftops for .eps or a working LaTeX/ghostscript install)
# filename = "homotopy_comparison.pdf" 
# savefig(final_plot, filename)
# println("\n=== Final Figure saved as $filename ===")

# Also display the plot to the screen for inspection
display(final_plot)