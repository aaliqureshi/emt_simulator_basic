module StaticInit

export run_static_init!

include("Models.jl"); using .Models
include("utils.jl"); using .Utils

function compute_shunt_currents(models)
    v_dp = @. models.bus.vd + 1im * models.bus.vq
    B_mat = build_B_matrix(models)
    i_cap = 1im * B_mat * v_dp
    i_cap_d = real(i_cap)
    i_cap_q = imag(i_cap)
    return i_cap_d, i_cap_q
end

function run_static_init!(sys)
    models = sys.models
    Y = sys.Y

    # compute bus power injections
    v_complex = @. models.bus.v * exp(1im * models.bus.theta)
    s_bus = v_complex .* conj(Y * v_complex)
    p_bus = real(s_bus)
    q_bus = imag(s_bus)

    # power balance: p_gen = p_load + p_network
    p_gen = zeros(length(models.bus.idx))
    p_gen[models.load.bus] += @. models.load.p
    p_gen += @. p_bus

    q_gen = zeros(length(models.bus.idx))
    q_gen[models.load.bus] += @. models.load.q
    q_gen += @. q_bus

    # compute line currents
    compute_line_currents!(models)

    # compute load currents
    load_id, load_iq = compute_load_currents(models)

    # compute shunt currents
    i_cap_d, i_cap_q = compute_shunt_currents(models)

    # compute generator and slack currents
    i_gen = zeros(Complex{Float64}, length(models.bus.idx))
    s_gen = @. p_gen[models.generator.bus] + 1im * q_gen[models.generator.bus]
    i_gen[models.generator.bus] = @. conj(s_gen / v_complex[models.generator.bus])

    s_slack = @. p_gen[models.slack.bus] + 1im * q_gen[models.slack.bus]
    i_slack = @. conj(s_slack / v_complex[models.slack.bus])

    i_slack_d = real(i_slack)
    i_slack_q = imag(i_slack)

    i_gen_d = real(i_gen)
    i_gen_q = imag(i_gen)

    # current balance verification
    incidence_matrix = sys.incidence_matrix
    i_d = incidence_matrix * models.line.i_d
    i_q = incidence_matrix * models.line.i_q

    i_d_balance = zeros(length(models.bus.idx))
    i_d_balance[models.generator.bus] += i_gen_d[models.generator.bus]
    i_d_balance[models.load.bus] -= load_id
    i_d_balance[:] -= i_cap_d
    i_d_balance[:] += i_d
    i_d_balance[models.slack.bus] += i_slack_d

    i_q_balance = zeros(length(models.bus.idx))
    i_q_balance[models.generator.bus] += i_gen_q[models.generator.bus]
    i_q_balance[models.load.bus] -= load_iq
    i_q_balance[:] -= i_cap_q
    i_q_balance[:] += i_q
    i_q_balance[models.slack.bus] += i_slack_q

    # compute generator internal voltage (E'q, delta)
    v_gen = @. (models.bus.vd[models.generator.bus] + 1im * models.bus.vq[models.generator.bus]) +
            1im*models.generator.x_d_prime*(i_gen_d[models.generator.bus] + 1im*i_gen_q[models.generator.bus])

    v_gen_mag = @. abs(v_gen)
    v_gen_angle = @. angle(v_gen)

    models.generator.e_q_prime[:] = @. v_gen_mag
    models.generator.delta[:] = @. v_gen_angle

    # transform generator currents to dq frame
    gen_id = @. i_gen_d[models.generator.bus] * sin(models.generator.delta) - i_gen_q[models.generator.bus] * cos(models.generator.delta)
    gen_iq = @. i_gen_d[models.generator.bus] * cos(models.generator.delta) + i_gen_q[models.generator.bus] * sin(models.generator.delta)

    models.generator.i_d[:] = @. gen_id
    models.generator.i_q[:] = @. gen_iq

    # set initial fault impedance (open circuit)
    models.fault.r_s = [1e10]
    models.fault.x_s = [1e10*2*pi*60]

    return nothing
end

end # module
