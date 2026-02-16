using Pkg
Pkg.activate(".")

include("../src/Models.jl")
using .Models

include("../src/data_loader.jl")
using .DataLoader

include("../src/utils.jl")
using .Utils


data_file = "cases/Simple_cases/test_cct.xlsx"

# function create_data()
#     bus1 = [1, 2, 4, 4, 4]
#     bus2 = [2, 3, 3, 2, 1]
#     gen_bus = [1, 4]
#     load_bus = [1, 3, 2]

#     gen_current = [5, 16]
#     load_current = [1, 20, 0]
# end

models = load_data(data_file) 


incidence_matrix = build_incidence_matrix(models)


models.load.bus
models.pv.bus
models.slack.bus
models.gencls.bus

id = zeros(length(models.bus.idx))
iq = zeros(length(models.bus.idx))



pv_current = [5]
slack_current = [16]
load_current = [1, 20 ,0]

gen_current = vcat(pv_current, slack_current)


id[models.gencls.bus] += @. gen_current
id[models.load.bus] -= @. load_current

line_current = [4, 13, 7, 9, 0]

id[:] += incidence_matrix * line_current