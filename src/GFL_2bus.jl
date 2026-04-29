using Pkg; Pkg.activate(".")
using MyDiffEq, Plots

# using OrdinaryDiffEq
include("MyTransformations.jl")
using .MyTransformations

function Voltage_Source(t)
    Vm = 480*sqrt(2)/sqrt(3)
    f = 60
    phi1 = 0
    phi2 = 2*pi/3
    phi3 = 4*pi/3
    V1 = Vm*cos(2*pi*f*t + phi1)
    V2 = Vm*cos(2*pi*f*t - phi2)
    V3 = Vm*cos(2*pi*f*t + phi2)
    return [V1, V2, V3]
end

# --- Parameters ---
const I_max = 300.0          # Current saturation limit (A, peak DQ magnitude)

const t_fault_on  = 30.0      # Fault applied at 3.0 s
const t_fault_off = 30.1      # Fault cleared at 3.1 s  (6 cycles at 60 Hz)
const R_fault = 0.09         # Fault shunt resistance at bus 2 (ohms)

function fx!(du, u, p, t)
    x1_pll, x2_pll, x1_ccd, x1_ccq, x1_pcd, x1_pcq, i_fa, i_fb, i_fc, i_12a, i_12b, i_12c, m_d, m_q, id_ref, iq_ref, v_2a, v_2b, v_2c = u

    wt = x2_pll

    # t = 0

    L_f = 100e-6
    if t < 10/60
        R_f = 1e13
    else
        R_f = 3e-3
    end

    R_12 = 1e-3
    L_12 = 100e-6

    f = 60
    Kp_pll = 0.785398
    Ki_pll = 592.176264

    Kp_pc = 0.0
    Ki_pc = 0.032139

    V1_abc = Voltage_Source(t)
    v_1a, v_1b, v_1c = V1_abc

    V2_abc = v_2a, v_2b, v_2c
    V2_dq0 = ABC_2_DQ0(V2_abc, wt)
    v_2d, v_2q, _ = V2_dq0

    M_dq0 = [m_d, m_q, 0]
    M_abc = DQ0_2_ABC(M_dq0, wt)
    m_a, m_b, m_c = M_abc

    m_a = clamp(m_a, -1.0, 1.0)
    m_b = clamp(m_b, -1.0, 1.0)
    m_c = clamp(m_c, -1.0, 1.0)

    I12_abc = [i_12a, i_12b, i_12c]
    If_abc = [i_fa, i_fb, i_fc]
    If_dq0 = ABC_2_DQ0(If_abc, wt)
    i_fd, i_fq, _ = If_dq0

    # --- FIX 1: Power reference schedule (proper elseif) ---
    if t < 15/60
        p_ref = 0.0
        q_ref = 0.0
    elseif t < 20/60
        p_ref = 250e3
        q_ref = 0.0
    elseif t < 2.5
        p_ref = 0.0
        q_ref = -10e1
    else
        # Pre-fault and during/post-fault: maintain last setpoint
        p_ref = 0.0
        q_ref = -10e1
    end

    Kp_cc = 0.02
    Ki_cc = 0.6

    V_dc = 800

    v_ta = m_a*V_dc/2
    v_tb = m_b*V_dc/2
    v_tc = m_c*V_dc/2

    v_td, v_tq, _ = ABC_2_DQ0([v_ta, v_tb, v_tc], wt)

    p_out = (v_2d*i_fd)*(3/2)
    q_out = (-v_2d*i_fq)*(3/2)

    # --- Current limiting (d-axis priority) ---
    id_ref_unsat = (p_ref - p_out)*Kp_pc + x1_pcd*Ki_pc
    iq_ref_unsat = (q_out - q_ref)*Kp_pc + x1_pcq*Ki_pc   # FIX 2: sign corrected

    # Clamp d-axis first, then q gets remaining margin
    id_ref_lim = clamp(id_ref_unsat, -I_max, I_max)
    iq_margin  = sqrt(max(I_max^2 - id_ref_lim^2, 0.0))
    iq_ref_lim = clamp(iq_ref_unsat, -iq_margin, iq_margin)

    # Detect saturation for anti-windup
    saturated = (id_ref_lim != id_ref_unsat) || (iq_ref_lim != iq_ref_unsat)

    # ---- Differential equations ----

    # PLL
    du[1] = v_2q
    du[2] = v_2q*Kp_pll + x1_pll*Ki_pll + 2π*f

    # Current controller integrators
    du[3] = id_ref - i_fd
    du[4] = iq_ref - i_fq

    # Power controller integrators (FIX 2: sign on du[6])
    if saturated
        # Anti-windup: freeze integrators when current is saturated
        du[5] = 0.0
        du[6] = 0.0
    else
        du[5] = p_ref - p_out
        du[6] = q_out - q_ref                  # FIX 2: flipped sign
    end

    # Filter inductor currents (abc)
    du[7]  = (v_ta - v_2a - R_f*i_fa) / L_f
    du[8]  = (v_tb - v_2b - R_f*i_fb) / L_f
    du[9]  = (v_tc - v_2c - R_f*i_fc) / L_f

    # Line currents (abc)
    du[10] = (v_1a - v_2a - R_12*i_12a) / L_12
    du[11] = (v_1b - v_2b - R_12*i_12b) / L_12
    du[12] = (v_1c - v_2c - R_12*i_12c) / L_12

    # ---- Algebraic equations ----

    # Modulation indices
    du[13] = ((id_ref - i_fd)*Kp_cc + x1_ccd*Ki_cc + v_2d - i_fq*2π*60*L_f)*(1/400) - m_d
    du[14] = ((iq_ref - i_fq)*Kp_cc + x1_ccq*Ki_cc + v_2q + i_fd*2π*60*L_f)*(1/400) - m_q

    # Current references (with saturation)
    du[15] = id_ref - id_ref_lim
    du[16] = iq_ref - iq_ref_lim

    # KCL at Bus 2 (with fault shunt)
    if t_fault_on <= t < t_fault_off
        i_fault_a = v_2a / R_fault
        i_fault_b = v_2b / R_fault
        i_fault_c = v_2c / R_fault
    else
        i_fault_a = 0.0
        i_fault_b = 0.0
        i_fault_c = 0.0
    end

    du[17] = i_12a + i_fa - i_fault_a
    du[18] = i_12b + i_fb - i_fault_b
    du[19] = i_12c + i_fc - i_fault_c
end

# --- Initial conditions ---
u0 = ones(19)
u0[17:19] = Voltage_Source(0)

# --- Mass matrix (1:12 differential, 13:19 algebraic) ---
M = zeros(19, 19)
for i in 1:12
    M[i, i] = 1
end

# --- Solve ---
dt = 5e-3
tspan = (0, 5.0)
tstops = []
# prob0 = ODEFunction(fx!, mass_matrix=M)
prob = ODEProblem(fx!, u0, tspan, (), M)
sol = Solve(prob, dt, method=:Euler, adaptive=false, tstops=tstops, always_new=true)



# ============================================================
#  Post-processing & Plots
# ============================================================
using Plots

wt   = [u[2]  for u in sol.u]
m_d  = [u[13] for u in sol.u]
m_q  = [u[14] for u in sol.u]
i_fa = [u[7]  for u in sol.u]
i_fb = [u[8]  for u in sol.u]
i_fc = [u[9]  for u in sol.u]

id_ref_sol = [u[15] for u in sol.u]
iq_ref_sol = [u[16] for u in sol.u]

v_2a = [u[17] for u in sol.u]
v_2b = [u[18] for u in sol.u]
v_2c = [u[19] for u in sol.u]

# DQ currents
i_sd = Float64[]
i_sq = Float64[]
for i in 1:length(sol.t)
    If_dq0 = ABC_2_DQ0([i_fa[i], i_fb[i], i_fc[i]], wt[i])
    push!(i_sd, If_dq0[1])
    push!(i_sq, If_dq0[2])
end

# Modulation in abc
m_a_sol = Float64[]
m_b_sol = Float64[]
m_c_sol = Float64[]
for i in 1:length(sol.t)
    m_abc = DQ0_2_ABC([m_d[i], m_q[i], 0.0], wt[i])
    push!(m_a_sol, clamp(m_abc[1], -1.0, 1.0))
    push!(m_b_sol, clamp(m_abc[2], -1.0, 1.0))
    push!(m_c_sol, clamp(m_abc[3], -1.0, 1.0))
end

# --- Current magnitude & saturation limit ---
i_mag = sqrt.(id_ref_sol.^2 .+ iq_ref_sol.^2)

p1 = plot(sol.t, i_mag, label="|I_ref|", xlabel="Time (s)", ylabel="Current (A)",
          title="Current Reference Magnitude vs Saturation Limit")
hline!(p1, [I_max], label="I_max = $I_max A", linestyle=:dash, color=:red)

# --- DQ currents ---
p2 = plot(sol.t, i_sd, label="i_d", xlabel="Time (s)", ylabel="Current (A)",
          title="DQ Currents")
plot!(p2, sol.t, i_sq, label="i_q")

# --- Bus 2 voltage (abc) ---
p3 = plot(sol.t, v_2a, label="v_2a", xlabel="Time (s)", ylabel="Voltage (V)",
          title="Bus 2 Voltage (abc)")
plot!(p3, sol.t, v_2b, label="v_2b")
plot!(p3, sol.t, v_2c, label="v_2c")

# --- Bus 2 voltage magnitude (DQ) ---
v_2d_sol = Float64[]
for i in 1:length(sol.t)
    V2_dq0 = ABC_2_DQ0((v_2a[i], v_2b[i], v_2c[i]), wt[i])
    push!(v_2d_sol, V2_dq0[1])
end
p4 = plot(sol.t, v_2d_sol, label="v_2d", xlabel="Time (s)", ylabel="Voltage (V)",
          title="Bus 2 DQ Voltage (v_2d ≈ magnitude)")

# --- Active & reactive power ---
p_out_sol = (3/2) .* v_2d_sol .* i_sd
q_out_sol = -(3/2) .* v_2d_sol .* i_sq

p5 = plot(sol.t, p_out_sol, label="P_out", xlabel="Time (s)", ylabel="Power",
          title="Active & Reactive Power")
plot!(p5, sol.t, q_out_sol, label="Q_out")

# --- Modulation indices ---
p6 = plot(sol.t, m_d, label="m_d", xlabel="Time (s)", ylabel="Modulation",
          title="DQ Modulation Indices")
plot!(p6, sol.t, m_q, label="m_q")

# --- Phase currents ---
p7 = plot(sol.t, i_fa, label="i_fa", xlabel="Time (s)", ylabel="Current (A)",
          title="Filter Currents (abc)")
plot!(p7, sol.t, i_fb, label="i_fb")
plot!(p7, sol.t, i_fc, label="i_fc")

# --- PLL angle ---
p8 = plot(sol.t, wt, label="θ_pll", xlabel="Time (s)", ylabel="Angle (rad)",
          title="PLL Angle")

display(plot(p1, p2, p3, p4, p5, p6, p7, p8, layout=(4, 2), size=(1200, 1600)))

# ============================================================
#  Small-Signal Analysis (Schur complement for DAE)
# ============================================================
using LinearAlgebra

function ss_analysis(t, u)
    n = length(u)
    n_diff = 12
    jac = zeros(n, n)
    ϵ = 1e-6
    du_eq = zeros(n)
    fx!(du_eq, u, nothing, t)
    u_pert = copy(u)
    for i in eachindex(u)
        u_pert[i] = u[i] + ϵ
        du_pert = zeros(n)
        fx!(du_pert, u_pert, nothing, t)
        jac[:, i] = (du_pert - du_eq) / ϵ
        u_pert[i] = u[i]
    end

    # Schur complement: reduce DAE to ODE eigenvalue problem
    f_x = jac[1:n_diff, 1:n_diff]
    f_y = jac[1:n_diff, (n_diff+1):end]
    g_x = jac[(n_diff+1):end, 1:n_diff]
    g_y = jac[(n_diff+1):end, (n_diff+1):end]

    As = f_x - f_y * inv(g_y) * g_x
    eig_vals = eigvals(As)
    return eig_vals
end

# Evaluate at a pre-fault steady-state point
t_eq = sol.t[end-10]
u_eq = sol.u[end-10]

eig_vals = ss_analysis(t_eq, u_eq)

p_eig = plot(real.(eig_vals), imag.(eig_vals), seriestype=:scatter,
             xlabel="Real", ylabel="Imaginary",
             title="Eigenvalues (Schur complement)", legend=false)
display(p_eig)
