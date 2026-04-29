using Pkg; Pkg.activate(".")
using MyDiffEq, Plots

function fx!(du, u, p, t)
    i, v_c, v_l = u

    v_m, r, c, l = p 

    ω = 2*pi*60

    v_in = v_m * sin(ω*t)

    # if t>=1.0
    #     s=1e-10
    # else
    #     s=1e10
    # end

    if 0.3 <= t <= 0.5
        s = 1e-10
    else
        s = 1e10
    end

    du[1] = v_l/l
    du[2] = i/c
    du[3] = v_in - i*r - i*s - v_l - v_c
end

u0 = [1.0, 1.0, 1.0]
v_m = 220e3
# r = 17.439
r = 0.0
c = 3.946e-6
l = 0.3e-3
p = (v_m, r, c, l)
dt = 1e-3

mass_matrix=zeros(3,3)
mass_matrix[1,1] = mass_matrix[2,2] = 1.0

prob = ODEProblem(fx!, u0, (0.0, 2.0), p, mass_matrix)
sol = Solve(prob, dt, method=:Euler, adaptive=false, tstops=[], always_new=true)

states = get_states(sol)

plot(sol.time, states[1,:])
plot(sol.time, states[2,:])
plot(sol.time, states[3,:])

t = 0:1e-3:0.2
vd = 0.7
vt = @. 15 * sin(2*pi*60*t)

i_s = 1

i = @. i_s * (exp(vd/vt) - 1)

plot(t, i)