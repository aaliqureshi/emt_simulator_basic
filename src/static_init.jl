# compute line currents
function compute_line_currents!(models)
    v_dp = @. models.bus.vd + 1im * models.bus.vq
    i_line = @. (v_dp[models.line.bus1_idx] - v_dp[models.line.bus2_idx]) / (models.line.R + 1im * models.line.X)
    models.line.i_d .= real(i_line)
    models.line.i_q .= imag(i_line)
end

compute_line_currents!(models)


function compute_load_currents(models)
    v_dp = @. models.bus.vd[models.load.bus] - 1im * models.bus.vq[models.load.bus]
    i_load = @. (models.load.p - 1im * models.load.q) / v_dp
    return real(i_load), imag(i_load)
end

load_id, load_iq = compute_load_currents(models)

# compute shunt currents
function compute_shunt_currents(models)
    v_dp = @. models.bus.vd + 1im * models.bus.vq
    B_mat = build_B_matrix(models)
    i_cap = 1im * B_mat * v_dp
    i_cap_d = real(i_cap)
    i_cap_q = imag(i_cap)
    return i_cap_d, i_cap_q
end

i_cap_d, i_cap_q = compute_shunt_currents(models)


i_gen = zeros(Complex{Float64}, length(models.bus.idx))
s_gen = @. p_gen[models.generator.bus] + 1im * q_gen[models.generator.bus]

i_gen[models.generator.bus] = @. conj(s_gen / v_complex[models.generator.bus])  

s_slack = @. p_gen[models.slack.bus] + 1im * q_gen[models.slack.bus]
i_slack = @. conj(s_slack / v_complex[models.slack.bus])

i_slack_d = real(i_slack)
i_slack_q = imag(i_slack)

i_gen_d = real(i_gen)
i_gen_q = imag(i_gen)

incidence_matrix = build_incidence_matrix(models)
i_d = incidence_matrix * models.line.i_d

i_q = incidence_matrix * models.line.i_q
# balance currents
i_d_balance = zeros(length(models.bus.idx))
i_d_balance[models.generator.bus] += i_gen_d[models.generator.bus]
i_d_balance[models.load.bus] -= load_id
i_d_balance[:] -= i_cap_d
i_d_balance[:] += i_d
i_d_balance[models.slack.bus] += i_slack_d


i_d_balance

# balance currents
i_q_balance = zeros(length(models.bus.idx))
i_q_balance[models.generator.bus] += i_gen_q[models.generator.bus]
i_q_balance[models.load.bus] -= load_iq
i_q_balance[:] -= i_cap_q
i_q_balance[:] += i_q
i_q_balance[models.slack.bus] += i_slack_q

i_q_balance


i_gen_d

# i_cap_d[models.generator.bus] - load_id

i_gen_d[models.generator.bus]

v_gen = @. (models.bus.vd[models.generator.bus] + 1im * models.bus.vq[models.generator.bus]) + 
        1im*models.generator.x_d_prime*(i_gen_d[models.generator.bus] + 1im*i_gen_q[models.generator.bus])

v_gen_mag = @. abs(v_gen)
v_gen_angle = @. angle(v_gen)

models.generator.e_q_prime[:] = @. v_gen_mag
models.generator.delta[:] = @. v_gen_angle

gen_id = @. i_gen_d[models.generator.bus] * sin(models.generator.delta) - i_gen_q[models.generator.bus] * cos(models.generator.delta)
gen_iq = @. i_gen_d[models.generator.bus] * cos(models.generator.delta) + i_gen_q[models.generator.bus] * sin(models.generator.delta)

models.generator.i_d[:] = @. gen_id
models.generator.i_q[:] = @. gen_iq

models.generator.x_d_prime