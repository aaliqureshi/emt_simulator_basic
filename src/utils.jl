module Utils

export build_Y_matrix, build_incidence_matrix, build_B_matrix

using Pkg
Pkg.activate(".")

#TODO: in future this function should update system object

function build_Y_matrix(models)
    # TODO: convert this to sparse matrices for better performance
    # TODO: add transformers
    n_bus = length(models.bus.idx)
    Y = zeros(ComplexF64, n_bus, n_bus)

    # create Y matrix
    for line_idx in eachindex(models.line.idx)
        Y[models.line.bus1_idx[line_idx], models.line.bus2_idx[line_idx]] = -1 / (models.line.R[line_idx] + 1im * models.line.X[line_idx])
        Y[models.line.bus2_idx[line_idx], models.line.bus1_idx[line_idx]] = -1 / (models.line.R[line_idx] + 1im * models.line.X[line_idx])
    end
    for row_idx in 1:n_bus
        Y[row_idx, row_idx] = Y[row_idx, row_idx] + sum(-1 .* Y[row_idx, :])
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

    num_bus = length(models.bus.idx)
    num_line = length(models.line.idx)
    incidence_matrix = zeros(Int32, num_bus, num_line)

    for line in collect(1:num_line)
        incidence_matrix[bus1_idx[line], line] = -1
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
    for (bus1, bus2, B) in zip(models.line.bus1_idx, models.line.bus2_idx, models.line.B)
        B_half = 1im * B / 2
        B_mat[bus1, bus1] += B_half
        B_mat[bus2, bus2] += B_half
    end

    # Add artficial capacitance to Every bus - this ensures mass matrix diagonal is not 0
    C_artifact = 1e-6
    B_artifact = 2*pi*60*C_artifact
    
    for i in eachindex(1:n_bus)
        B_mat[i, i] += 1im * B_artifact
    end

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