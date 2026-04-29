module MyTransformations

export ABC_2_αβ0, ABC_2_DQ0, DQ0_2_αβ0, DQ0_2_ABC

function ABC_2_αβ0(V1, V2, V3)
    T = [1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2; 1/2 1/2 1/2]
    Vα, Vβ, V0 = (2/3)*T*[V1, V2, V3]
    return [Vα, Vβ, V0]
end

function ABC_2_DQ0(Vabc, wt)
    Vα, Vβ, V0 = ABC_2_αβ0(Vabc...)
    theta = wt

    ## aligned with A
    # Vd = Vα*cos(theta) + Vβ*sin(theta)
    # Vq = -Vα*sin(theta) + Vβ*cos(theta)

    ## 90 deg behind A
    Vd = Vα*sin(wt) - Vβ*cos(wt)
    Vq = Vα*cos(wt) + Vβ*sin(wt)

    V0 = V0
    return [Vd, Vq, V0]
end

function DQ0_2_αβ0(Vdq0, θ)
    Vd, Vq, V0 = Vdq0

    ## aligned with A
    # T = [cos(θ) -sin(θ); sin(θ) cos(θ)]

    ## 90 deg behind A
    T = [sin(θ) cos(θ); -cos(θ) sin(θ)]


    Vα, Vβ = T*[Vd, Vq]
    V0 = V0
    return [Vα, Vβ, V0]
end

function αβ0_2_ABC(Vα, Vβ, V0)
    T = [1 0 1; -1/2 sqrt(3)/2 1; -1/2 -sqrt(3)/2 1]
    V1, V2, V3 = T*[Vα, Vβ, V0]
    return [V1, V2, V3]
end

function DQ0_2_ABC(Vdq0, θ)
    Vα, Vβ, V0 = DQ0_2_αβ0(Vdq0, θ)
    Vabc = αβ0_2_ABC(Vα, Vβ, V0)
    return Vabc
end


end # module







