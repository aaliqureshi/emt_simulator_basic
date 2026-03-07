module Utils

export build_Y_matrix, build_incidence_matrix, build_B_matrix

using Pkg
Pkg.activate(".")

#TODO: in future this function should update system object

function build_Y_matrix(models)
    # TODO: convert this to sparse matrices for better performance
    n_bus = length(models.bus.idx)
    Y = zeros(ComplexF64, n_bus, n_bus)

    # create Y matrix (with transformer tap handling)
    for line_idx in eachindex(models.line.idx)
        i = models.line.bus1_idx[line_idx]
        j = models.line.bus2_idx[line_idx]
        y = 1 / (models.line.R[line_idx] + 1im * models.line.X[line_idx])
        t = models.line.tap[line_idx]

        Y[i, j] += -y / t
        Y[j, i] += -y / t
        Y[i, i] += y / (t^2)
        Y[j, j] += y
    end

    B_mat = _build_B_matrix(models)

    # update Y matrix
    Y = Y + B_mat

    return Y
end


## add incidence matrix here
function build_incidence_matrix(models)
    bus1_idx = models.line.bus1_idx
    bus2_idx = models.line.bus2_idx
    tap = models.line.tap

    num_bus = length(models.bus.idx)
    num_line = length(models.line.idx)
    # transformer convention: Bus i ---(a:1)---[Z]--- Bus j
    # current at bus i is i_line/a, current at bus j is i_line
    incidence_matrix = zeros(Float64, num_bus, num_line)

    for line in collect(1:num_line)
        incidence_matrix[bus1_idx[line], line] = -1 / tap[line]
        incidence_matrix[bus2_idx[line], line] = 1
    end

    return incidence_matrix
end


function _build_B_matrix(models)
    """
    Build B matrix.
    """
    n_bus = length(models.bus.idx)
    B_mat = zeros(ComplexF64, n_bus, n_bus)
    # transformer convention: Bus i ---(a:1)---[B/2]---[Z]---[B/2]--- Bus j
    # bus i shunt sees V_i/a, so admittance referred to bus i is B/(2a²)
    for (bus1, bus2, B, a) in zip(models.line.bus1_idx, models.line.bus2_idx, models.line.B, models.line.tap)
        B_mat[bus1, bus1] += 1im * B / (2 * a^2)
        B_mat[bus2, bus2] += 1im * B / 2
    end

    # Add artficial capacitance to Every bus - this ensures mass matrix diagonal is not 0
    # C_artifact = 1e-6
    # B_artifact = 2*pi*60*C_artifact
    
    # for i in eachindex(1:n_bus)
    #     B_mat[i, i] += 1im * B_artifact
    # end

    return B_mat
end


function build_B_matrix(models)
    """
    This function returns the absolute value of the B matrix.
    """
    B_mat = _build_B_matrix(models)
    return abs.(B_mat)
end


end #module