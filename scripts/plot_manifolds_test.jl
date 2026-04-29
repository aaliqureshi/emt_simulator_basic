# using Plots

# g1(x) = x^2
# g2(x) = x^2 + 2x

# g3(x,y) = x^2 + y^2 - 1

# x = -5:0.1:5

# y = -5:0.1:5
# m1 = @. g1(x)
# m2 = @. g2(x)

# m3 = @. g3(x,y)

# plot(x, m1)
# plot!(x, m2)

# plot3d(x,y,m3)


# # Define a grid of x and y values
# x = range(-1.5, 1.5, length=100)
# y = range(-1.5, 1.5, length=100)

# # Define the function z = x^2 + y^2
# z(x, y) = x^2 + y^2

# # Plot the contour where z = 1
# contour(x, y, z, levels=[1], 
#         aspect_ratio=:equal, 
#         colorbar=false, 
#         title="Contour Plot: x² + y² = 1")

using OrdinaryDiffEq
using Plots

r = 10
function fx!(du,u,p,t)
    # if t <= 1.0
    #     r = 32
    # else
    #     r = 10
    # end
    x, y = u
    # r = p[1]
    du[1] = -x + cos(y)
    du[2] = x^2 + y^2 - r
end

mass_matrix = zeros(2,2)
mass_matrix[1,1] = 1.0
u0 = [3,1]
u0 = [4,-4]
# u0=[0.718,5.61]
prob0 = ODEFunction(fx!, mass_matrix = mass_matrix)
prob = OrdinaryDiffEq.ODEProblem(prob0, u0, [0, 2.0], ())
sol = OrdinaryDiffEq.solve(prob, ImplicitEuler(;nlsolve=NLNewton(always_new=true)), adaptive=false, dt=1e-3)

X = [u[1] for u in sol.u]
Y = [u[2] for u in sol.u]

# Plot the algebraic manifold: x² + y² = 1 (unit circle)
θ = range(0, 2π, length=200)

# r=32
x_circle = sqrt(r).*cos.(θ)
y_circle = sqrt(r).*sin.(θ)

p1 = plot(x_circle, y_circle, label="Manifold: x² + y² = $r",
     linewidth=2, linestyle=:dash, color=:blue,
     aspect_ratio=:equal, xlabel="x", ylabel="y",
     title="DAE Solution on Algebraic Manifold")

# Overlay the DAE solution trajectory
plot!(p1, X, Y, label="DAE solution (x(t), y(t))",
     linewidth=2, color=:red)
scatter!(p1, [u0[1]], [u0[2]], label="Initial point", markersize=6, color=:green)
scatter!(p1, [X[end]], [Y[end]], label="Final point", markersize=6, color=:orange)

# Plot the algebraic residual g(x,y) = x² + y² - 1 over time
g_residual = X.^2 .+ Y.^2 .- r
p2 = plot(sol.t, g_residual, label="g(x,y) = x² + y² - 1",
     linewidth=2, xlabel="t", ylabel="Algebraic residual",
     title="Algebraic Constraint Residual Over Time")
hline!(p2, [0.0], linestyle=:dash, color=:black, label="")

p = plot(p1, p2, layout=(1,2), size=(1000, 450))
display(p)
