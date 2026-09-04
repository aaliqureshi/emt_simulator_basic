module BusModel

export Bus, balance!, phasor2DP!

mutable struct Bus{T<:Real}
    idx :: Vector{Int32}
    orig_idx :: Vector{Int32}
    v :: Vector{T}
    theta :: Vector{T}
    vd:: Vector{T}
    vq:: Vector{T}

    function Bus{T}(n::Integer) where {T<:Real}
        new{T}(Vector{Int32}(undef, n),
               Vector{Int32}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n),
               Vector{T}(undef, n))
    end
end

function phasor2DP!(bus)
    vdq = @. bus.v * exp(1im * bus.theta)
    bus.vd = real(vdq)
    bus.vq = imag(vdq)
end

function balance!(du, u, p)
    T = eltype(u)
    address, models, incidence_matrix, C_eq, non_slack_buses, _ = p

    bus = models.bus
    line = models.line
    generator = models.generator
    fault = models.fault
    load = models.load
    slack = models.slack

    bus_vd = Vector{T}(bus.vd)
    bus_vq = Vector{T}(bus.vq)

    gen_delta = u[address["delta"]]
    gen_omega = u[address["omega"]]
    gen_id = u[address["gen_id"]]
    gen_iq = u[address["gen_iq"]]
    line_id = u[address["line_id"]]
    line_iq = u[address["line_iq"]]

    bus_vd[non_slack_buses] = @. u[address["balance_d"]]
    bus_vq[non_slack_buses] = @. u[address["balance_q"]]

    fault_id = u[address["fault_id"]]
    fault_iq = u[address["fault_iq"]]

    id = zeros(T, length(bus.idx))
    iq = zeros(T, length(bus.idx))

    # ZIP load model — computed directly in (d,q) to avoid ForwardDiff singularities.
    #
    # Key identity: conj(S/V) = (P·vd + Q·vq + j(P·vq - Q·vd)) / |V|²
    # This lets us express load current using |V|² = vd²+vq² (polynomial, no sqrt)
    # instead of S/V (which has a 1/V singularity ForwardDiff can't cancel).
    #
    # Z component: I = conj(S₀)·V / |V₀|²     — linear in V, no singularity
    # I component: I = conj(S₀)·V / (|V₀|·|V|) — regularize |V| at low voltage
    # P component: I = conj(S₀)·V / |V|²       — regularize |V|² at low voltage

    # k_pz = T(0.333)
    # k_pi = T(0.333)
    # k_pp = T(0.333)

    k_pz = T(0.7)
    k_pi = T(0.1)
    k_pp = T(0.2)
    v_min = T(0.7)   # voltage floor (pu) — I and P components freeze below this

    vd_l = bus_vd[load.bus]
    vq_l = bus_vq[load.bus]
    p_l  = Vector{T}(load.p)
    q_l  = Vector{T}(load.q)

    # |V₀|² from initial operating point (constant w.r.t. ForwardDiff)
    v0_sq = @. bus.vd[load.bus]^2 + bus.vq[load.bus]^2

    # |V|² from current state — polynomial in (vd, vq), ForwardDiff-smooth everywhere
    v_sq = @. vd_l^2 + vq_l^2

    # ── Z component: I_z = conj(S₀) · V / |V₀|²  (exact Ohm's law, no singularity) ──
    i_d_z = @. k_pz * (p_l * vd_l + q_l * vq_l) / v0_sq
    i_q_z = @. k_pz * (p_l * vq_l - q_l * vd_l) / v0_sq

    # ── I component: I_i = conj(S₀) · V / (|V₀| · |V|) ──
    # Regularize: replace |V|² with max(|V|², v_min²), then |V| = sqrt(that)
    v_sq_i = @. max(v_sq, v_min^2)
    # Denominator: |V₀| · |V_eff| = sqrt(|V₀|² · v_sq_eff)  — one sqrt of a clamped value
    denom_i = @. sqrt(v0_sq * v_sq_i)
    i_d_i = @. k_pi * (p_l * vd_l + q_l * vq_l) / denom_i
    i_q_i = @. k_pi * (p_l * vq_l - q_l * vd_l) / denom_i

    # ── P component: I_p = conj(S₀) · V / |V|² ──
    # Regularize: replace |V|² with max(|V|², v_min²)
    v_sq_p = @. max(v_sq, v_min^2)
    i_d_p = @. k_pp * (p_l * vd_l + q_l * vq_l) / v_sq_p
    i_q_p = @. k_pp * (p_l * vq_l - q_l * vd_l) / v_sq_p

    # Total load current
    i_load_d = @. i_d_z + i_d_i + i_d_p
    i_load_q = @. i_q_z + i_q_i + i_q_p

    id[generator.bus] += @. gen_id * sin(gen_delta) + gen_iq * cos(gen_delta)
    id[load.bus] -= i_load_d
    id[fault.bus] -= @. fault_id
    id[:] += incidence_matrix * line_id

    iq[generator.bus] += @. gen_iq * sin(gen_delta) - gen_id * cos(gen_delta)
    iq[load.bus] -= i_load_q
    iq[fault.bus] -= @. fault_iq
    iq[:] += incidence_matrix * line_iq

    omega = 2*pi*60
    w1 = @. omega*C_eq * bus_vq
    w2 = @. omega*C_eq * bus_vd
    id[:] += w1
    iq[:] -= w2

    p_h = zeros(T, length(bus.idx))
    q_h = zeros(T, length(bus.idx))

    p_h[non_slack_buses] = @. id[non_slack_buses] * bus_vd[non_slack_buses] + iq[non_slack_buses] * bus_vq[non_slack_buses]
    q_h[non_slack_buses] = @. id[non_slack_buses] * bus_vq[non_slack_buses] - iq[non_slack_buses] * bus_vd[non_slack_buses]

    # du[address["balance_d"]] = @. id[non_slack_buses]
    # du[address["balance_q"]] = @. iq[non_slack_buses]

    du[address["balance_d"]] = @. p_h[non_slack_buses]
    du[address["balance_q"]] = @. q_h[non_slack_buses]
end

end # module
