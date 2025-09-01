module Utils

export build_Y_matrix, build_incidence_matrix

using Pkg
Pkg.activate(".")

#TODO: in future this function should update system

function build_Y_matrix(models)
    n_bus = length(models.bus.idx)
    Y = zeros(ComplexF64, n_bus, n_bus)

    # create Y matrix
    #TODO: add the B matrix
    for line_idx in eachindex(models.line.idx)
        Y[models.line.bus1_idx[line_idx], models.line.bus2_idx[line_idx]] = -1 / (models.line.R[line_idx] + 1im * models.line.X[line_idx])
        Y[models.line.bus2_idx[line_idx], models.line.bus1_idx[line_idx]] = -1 / (models.line.R[line_idx] + 1im * models.line.X[line_idx])
    end
    for row_idx in 1:n_bus
        Y[row_idx, row_idx] = Y[row_idx, row_idx] + sum(-1 .* Y[row_idx, :])
    end

    return abs.(Y), angle.(Y)
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


end #module