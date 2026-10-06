using LinearAlgebra, Plots


# const V_ref  = 1.0
# const R_g    = 0.1
# const P_load = 0.9
# const I_max  = 3.0
# const V_init = 0.9

# const V_ref  = 1.0
# const R_g    = 0.1
# const P_load = 0.9
# const I_max  = 1.5
# const V_init = 0.9
# const I_load = 0.7
# const r_fault = 0.01
# const kp = 2.0
# const x = 0.0

const I_max = 1.0
const I_load = 0.8
const kp = 1.0
const x = 0.0
const V_ref = 1.0
V_init = 0.99


# function fx(y,μ)
#     return y^2 - y + μ/4
# end

# function dfx(y)
#     return 2*y - 1
# end

function fx(y,μ)
    # return I_max - P_load / y
    # return y*I_max - P_load
    # return sqrt(I_max^2 - (x+kp*(V_ref-y))^2) - I_load
    # return sqrt(I_max^2 - (min(I_max, abs(x+kp*(V_ref-y))))^2) - I_load
    # return sqrt(I_max^2 - clamp(x+kp*(V_ref-y), -I_max, I_max)^2) - y/r_fault
    # return I_max^2 - I_load^2 - (x+kp*(V_ref - y))^2
    # return (y-0.4)*(y-500)
    # nm = I_max
    # dn = sqrt(1 + ((x+kp*(V_ref-y))/I_max)^2)
    # return nm/dn - I_load - y/r_fault
    i_q = clamp(x + kp*(V_ref - y), -I_max, I_max)
    return sqrt(I_max^2 - i_q^2) - I_load
end

function dfx(y)
    # return P_load / y^2 
    # return I_max
    # return -1/R_g 
    # return (kp^2*(V_ref-y))/(sqrt(I_max^2 - kp^2*(V_ref-y)^2)) - 1/r_fault
    # return 2*kp*(x + kp*(V_ref-y))
    # return 2*kp^2*(V_ref-y)
    # nm = kp*(x+kp*(V_ref-y))
    # dn = I_max * (1 + ((x+kp*(V_ref-y))/(I_max))^2)^(3/2)
    # dn = sqrt(I_max^2 - (x + kp*(V_ref - y))^2)
    # return nm/dn - (1/r_fault)
    # return nm/dn
    if -I_max < x + kp*(V_ref - y) < I_max
        nm = kp*(x + kp*(V_ref-y))
        dn = sqrt(I_max^2 - (x + kp*(V_ref - y))^2)
        return nm/dn
    else
        return 0
    end
    return nm/dn
end

function NR(x_k, μ)
    return x_k - fx(x_k, μ)/(dfx(x_k))
end

function HBNR(x_k, x_prev, μ)
    beta = 0.1
    return x_k*(1+beta) - beta*x_prev - fx(x_k, μ)/(dfx(x_k))
end

function DNR(x_k, μ)
    alpha = 0.1
    return x_k - alpha * fx(x_k, μ)/(dfx(x_k))
end

function LM(x_k, μ)
    return x_k - fx(x_k, μ)*(dfx(x_k))
end

# function Reg_LM(x_k, μ)
#     fac = 0.01
#     return x_k - fx(x_k, μ)*(dfx(x_k) + fac)
# end

function Reg_LM(x_k, μ)
    alpha = 0.01
    fac = 5.0
    # return x_k - alpha * fx(x_k, μ)
    # return x_k - alpha * fx(x_k, μ) * dfx(x_k)
    return x_k - alpha * fx(x_k, μ) * (dfx(x_k) + fac)
end


function main(solver)
    iter = 0
    max_iter = 1000
    converged = false
    # guess_prev = guess = 0.9
    guess_prev = guess = V_init
    μ = 1.0
    # solver=:nr
    res_fx = []
    res_dfx = []
    guess_hist=[]
    push!(res_fx, fx(guess_prev, μ))
    push!(res_dfx, dfx(guess_prev))
    push!(guess_hist, guess)

    while iter<max_iter && !converged
        println("iter: $iter, residual: $(res_fx[end]), guess: $(guess)")
        if solver === :nr
            guess = NR(guess_prev, μ)
        elseif solver ===:lm
            guess = LM(guess_prev, μ)
        elseif solver ===:rlm
            guess = Reg_LM(guess_prev, μ)
        elseif solver ===:dnr
            guess = DNR(guess_prev, μ)
        elseif solver ===:hbnr
            if iter == 0
                x_prev = 0.0
            else
                x_prev = guess_hist[iter]
            end
            guess = HBNR(guess_prev, x_prev, μ)
        end
        push!(res_fx, fx(guess, μ))
        push!(res_dfx, dfx(guess))
        push!(guess_hist, guess)
        # @show abs(guess_prev - guess)
        guess_prev = guess

        if abs(res_fx[end]) < 1e-7
            converged = true
            println("converged")
        end

        iter+=1

    end
    return res_fx, res_dfx, guess_hist
end

solver=:nr
# solver = :lm
# solver = :rlm
# solver=:dnr
# solver=:hbnr
res_fx, res_dfx, guess = main(solver)

using Plots

plot(res_fx)
plot(res_dfx)
plot!(diff(guess))
plot(guess)
guess



# const V_ref  = 1.0
# const R_g    = 0.1
# const P_load = 0.9
# const I_max  = 1.0
# const V_init = 0.8
# const I_load = 0.9
# const r_fault = 0.005
# const kp = 1.0
# const x = 0.001

j = -2.0:0.1:2.0
# j = 0.1:0.01:1.0
vv = @. fx(j, 1.0)
kk = @. dfx(j)


plot(j, vv)
plot(j, kk)

dfx(0.8)

NR(1.02, 1.0)
NR(1.43, 1.0)

fx(0.91)