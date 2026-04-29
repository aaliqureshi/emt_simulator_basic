module SteadyStateAnalysis

using ForwardDiff
using LinearAlgebra
using Printf
using Dates

export small_signal_analysis, write_small_signal_report

# Build a pure residual f(u) -> du from the in-place DAE function, with t = 0
# and a captured parameter tuple. The fault is expected to be off at t = 0
# (caller is responsible for pushing fault times outside the analysis window).
function _build_residual(solve_func, p)
    function f(u::AbstractVector{T}) where {T}
        du = zeros(T, length(u))
        solve_func(du, u, p, zero(T))
        return du
    end
    return f
end

# Human-readable label per state index, derived from the address dictionary.
function _state_labels(address, models, non_slack_buses)
    nstates = maximum(maximum.(values(address)))
    labels = fill("", nstates)

    gen_bus = models.generator.bus
    for (i, k) in enumerate(address["delta"]);   labels[k] = "delta_gen_bus$(gen_bus[i])";  end
    for (i, k) in enumerate(address["omega"]);   labels[k] = "omega_gen_bus$(gen_bus[i])";  end
    for (i, k) in enumerate(address["gen_id"]);  labels[k] = "i_d_gen_bus$(gen_bus[i])";    end
    for (i, k) in enumerate(address["gen_iq"]);  labels[k] = "i_q_gen_bus$(gen_bus[i])";    end

    line = models.line
    for (i, k) in enumerate(address["line_id"])
        labels[k] = "i_d_line$(line.idx[i])(b$(line.bus1_idx[i])->b$(line.bus2_idx[i]))"
    end
    for (i, k) in enumerate(address["line_iq"])
        labels[k] = "i_q_line$(line.idx[i])(b$(line.bus1_idx[i])->b$(line.bus2_idx[i]))"
    end

    fault_bus = models.fault.bus
    for (i, k) in enumerate(address["fault_id"]); labels[k] = "i_d_fault_bus$(fault_bus[i])"; end
    for (i, k) in enumerate(address["fault_iq"]); labels[k] = "i_q_fault_bus$(fault_bus[i])"; end

    for (i, k) in enumerate(address["balance_d"]); labels[k] = "vd_bus$(non_slack_buses[i])"; end
    for (i, k) in enumerate(address["balance_q"]); labels[k] = "vq_bus$(non_slack_buses[i])"; end

    return labels
end

# Pair left eigenvectors to right eigenvectors by closest matching eigenvalue.
function _match_left_right(lam_R, lam_L, V_L)
    n = length(lam_R)
    W = similar(V_L, eltype(V_L), size(V_L))
    used = falses(n)
    for i in 1:n
        if !isfinite(lam_R[i])
            W[:, i] .= NaN
            continue
        end
        best_j, best_d = 0, Inf
        for j in 1:n
            used[j] && continue
            isfinite(lam_L[j]) || continue
            # d = abs(lam_L[j] - lam_R[i])
            d = abs(lam_L[j] - conj(lam_R[i]))
            if d < best_d
                best_d = d
                best_j = j
            end
        end
        if best_j > 0
            used[best_j] = true
            W[:, i] = V_L[:, best_j]
        else
            W[:, i] .= NaN
        end
    end
    return W
end

"""
    small_signal_analysis(solve_func, sys, address, mass_matrix, u0;
                          case_file=nothing, report_dir=".", ss_tol=1e-6, top_n=5)

Linearize the mass-matrix DAE `M ⋅ du/dt = solve_func(du, u, p, t=0)` about `u0`,
solve the generalized eigenproblem `J v = λ M v`, and return a NamedTuple with
eigenvalues, right/left eigenvectors, participation factors, and the residual
norm at `u0`. If `case_file` is given, also write a `.txt` report named after
the case file into `report_dir`.

The fault model is forced off during linearization by pushing `t_fault`/`t_clear`
to `Inf`; the original values are restored on exit.
"""
function small_signal_analysis(solve_func, sys, address, mass_matrix, u0;
                               case_file::Union{Nothing,AbstractString}=nothing,
                               report_dir::AbstractString=".",
                               ss_tol::Real=1e-6,
                               top_n::Integer=5)
    saved_t_fault = copy(sys.models.fault.t_fault)
    saved_t_clear = copy(sys.models.fault.t_clear)
    sys.models.fault.t_fault .= Inf
    sys.models.fault.t_clear .= Inf

    try
        lambda = 1.0
        p = (address, sys.models, sys.incidence_matrix, sys.C_eq, sys.non_slack_buses, lambda)

        f = _build_residual(solve_func, p)

        # 1) verify u0 is at steady state
        r0 = f(u0)
        residual_norm = norm(r0)
        is_ss = residual_norm <= ss_tol
        if !is_ss
            @warn "Initial condition is not at steady state — proceeding anyway" residual_norm tol=ss_tol
        end
        u_ss = copy(u0)

        # 2) Jacobian
        J = ForwardDiff.jacobian(f, u_ss)

        # 3) Generalized eigenproblem  J v = λ M v
        F_R = eigen(J, mass_matrix)
        lam = F_R.values
        V   = F_R.vectors

        # left eigenvectors satisfy  Jᵀ w = λ Mᵀ w
        F_L = eigen(Matrix(J'), Matrix(mass_matrix'))
        W   = _match_left_right(lam, F_L.values, F_L.vectors)

        # @show (W * N) 

        # 4) Participation factors  p_ki = w_i[k] v_i[k] / (w_iᵀ M v_i)
        nstates = size(J, 1)
        n = length(lam)
        P = fill(complex(NaN), nstates, n)
        for i in 1:n
            isfinite(lam[i]) || continue
            w_i = W[:, i]
            v_i = V[:, i]
            # denom = (W[:, i]' * mass_matrix * V[:, i])
            denom = w_i' * mass_matrix * v_i
            # denom = sum(abs.(w_i) .* abs.(v_i))
            row_weight = w_i' * mass_matrix
            if abs(denom) < 1e-12
                continue
            end

        #     abs(denom) < 1e-14 && continue
            for k in 1:nstates
                # P[k, i] = (abs(W[k, i]) * abs(V[k, i])) / denom
                P[k, i] = v_i[k] * row_weight[k] / denom
            
            end
        end

        # for i in 1:nstates
        #     for j in 1:n
        #     isfinite(lam[j]) || continue
        #     aa = abs(W[i, j]) * abs(V[j, i])
        #     denom = @. sum(abs(W[j, 1:2]) * abs(V[1:2,j]))
        #     P[i, j] = aa / denom
        #     end
        # end

        

        labels = _state_labels(address, sys.models, sys.non_slack_buses)

        result = (
            eigenvalues   = lam,
            right_vectors = V,
            left_vectors  = W,
            participation = P,
            u_ss          = u_ss,
            residual_norm = residual_norm,
            is_steady_state = is_ss,
            jacobian      = J,
            mass_matrix   = mass_matrix,
            state_labels  = labels,
        )

        if case_file !== nothing
            base = splitext(basename(case_file))[1]
            isdir(report_dir) || mkpath(report_dir)
            report_path = joinpath(report_dir, base * "_small_signal.txt")
            write_small_signal_report(report_path, result; top_n=top_n, case_file=case_file)
            @info "Small-signal report written" report_path
        end

        return result
    finally
        sys.models.fault.t_fault .= saved_t_fault
        sys.models.fault.t_clear .= saved_t_clear
    end
end

"""
    write_small_signal_report(path, result; top_n=5, case_file="")

Write a human-readable summary of `small_signal_analysis` output to `path`:
header, eigenvalue table sorted by real part (least damped first), algebraic
(infinite) modes flagged, and top-`top_n` participating states per finite mode.
"""
function write_small_signal_report(path::AbstractString, result; top_n::Integer=5, case_file::AbstractString="")
    lam     = result.eigenvalues
    labels  = result.state_labels
    P       = result.participation
    finite  = isfinite.(lam)
    n_fin   = count(finite)
    n_inf   = count(.!finite)

    idx_finite = findall(finite)
    sort!(idx_finite, by = i -> -real(lam[i]))

    open(path, "w") do io
        println(io, "Small-Signal Analysis Report")
        println(io, repeat("=", 72))
        println(io, "Case file        : ", case_file)
        println(io, "Generated        : ", string(now()))
        println(io, "States           : ", length(labels))
        println(io, "Linearization t  : 0.0  (fault off)")
        @printf(io, "Steady-state |f| : %.3e  (%s)\n",
                result.residual_norm,
                result.is_steady_state ? "OK" : "WARN: not at steady state")
        println(io)

        println(io, "Eigenvalue summary")
        println(io, repeat("-", 72))
        @printf(io, "  finite eigenvalues    : %d\n", n_fin)
        @printf(io, "  algebraic (Inf) modes : %d\n", n_inf)
        println(io)

        println(io, "Finite eigenvalues  (sorted by real part, least damped first)")
        println(io, repeat("-", 72))
        @printf(io, "%5s %18s %18s %12s %12s\n",
                "#", "real", "imag", "freq[Hz]", "damping")
        for (rank, i) in enumerate(idx_finite)
            re  = real(lam[i])
            im_ = imag(lam[i])
            freq = abs(im_) / (2*pi)
            mag  = sqrt(re^2 + im_^2)
            zeta = mag > 0 ? -re/mag : 0.0
            @printf(io, "%5d %18.6e %18.6e %12.4f %12.4f\n",
                    rank, re, im_, freq, zeta)
        end
        println(io)

        if n_inf > 0
            println(io, "Algebraic (infinite) modes")
            println(io, repeat("-", 72))
            for (k, i) in enumerate(findall(.!finite))
                @printf(io, "  alg %4d : λ = %s\n", k, string(lam[i]))
            end
            println(io)
        end

        println(io, "Participation factors  (top $(top_n) per finite mode, |p_ki| normalized)")
        println(io, repeat("-", 72))
        for (rank, i) in enumerate(idx_finite)
            re  = real(lam[i])
            im_ = imag(lam[i])
            p_i = abs.(P[:, i])
            s   = sum(p_i)
            p_n = s > 0 ? p_i ./ s : p_i
            order = sortperm(p_n, rev=true)
            @printf(io, "Mode %4d   λ = %.6e %+.6ej\n", rank, re, im_)
            ntop = min(top_n, length(order))
            for k in 1:ntop
                idx = order[k]
                @printf(io, "    %-44s  p = %.4f\n", labels[idx], p_n[idx])
            end
            println(io)
        end
    end
    return path
end

end # module
