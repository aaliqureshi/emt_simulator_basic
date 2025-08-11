module DataLoader

export load_data, models_dict

using Pkg
Pkg.activate(".")

using XLSX
using DataFrames

include("Models.jl")
using .Models


models_dict = Dict{String, Any}()

function load_data(data_file::String)
    if isfile(data_file)
        xf = XLSX.readxlsx(data_file)

        for sheet in XLSX.sheetnames(xf)
            df = DataFrame(XLSX.readtable(data_file, sheet))
            if sheet == "Bus"
                bus = Bus{Float64}(length(df.idx))
                bus.idx = Int32.(df.idx)
                bus.v = Float64.(df.v0)
                bus.theta = Float64.(df.a0)
                bus.vd = zeros(Float64, length(df.idx))
                bus.vq = zeros(Float64, length(df.idx))
                bus.i_d = zeros(Float64, length(df.idx))
                bus.i_q = zeros(Float64, length(df.idx))
                models_dict["bus"] = bus
            end
            if sheet == "Line"
                line = Line{Float64}(length(df.idx))
                line.idx = Int32.(df.idx)
                line.bus1_idx = Int32.(df.bus1)
                line.bus2_idx = Int32.(df.bus2)
                line.R = Float64.(df.r)
                line.X = Float64.(df.x)
                line.i_d = zeros(Float64, length(df.idx))
                line.i_q = zeros(Float64, length(df.idx))
                models_dict["line"] = line
            end
            if sheet == "GENCLS"
                generator = Generator{Float64}(length(df.idx))
                generator.idx = Int32.(df.idx)
                generator.bus = Int32.(df.bus)
                generator.delta = ones(Float64, length(df.idx))
                generator.omega = ones(Float64, length(df.idx))
                generator.M = Float64.(df.M)
                # generator.p_m = Float64.(df.p_m)
                generator.i_d = zeros(Float64, length(df.idx))
                generator.i_q = zeros(Float64, length(df.idx))
                generator.x_d_prime = Float64.(df.xd1)
                generator.e_q_prime = ones(Float64, length(df.idx))
        
                df_pv = DataFrame(XLSX.readtable(data_file, "PV"))
                generator.p_m = Float64.(df_pv.p0)
                generator.q_m = Float64.(df_pv.q0)
        
                models_dict["generator"] = generator
            end
        
            if sheet == "Fault"
                fault = Fault{Float64}(length(df.idx))
                fault.bus = Int32.(df.bus)
                fault.r_s = Float64.(df.rf)
                fault.l_s = Float64.(df.xf)
        
                models_dict["fault"] = fault
            end
            if sheet == "PQ"
                load = Load{Float64}(length(df.idx))
                load.bus = Int32.(df.bus)
                load.p = Float64.(df.p0)
                load.q = Float64.(df.q0)
        
                models_dict["load"] = load
            end
        
            if sheet == "Slack"
                slack = Slack(length(df.idx))
                slack.bus = Int32.(df.bus)
        
                models_dict["slack"] = slack
            end
        end

        return models_dict
    else
        error("File $data_file does not exist")
    end
end

# data_file = "cases/SMIB_Chow/SMIB_RL_Line_DrCui.xlsx"

# load_data(data_file)

# xf = XLSX.readxlsx(data_file)

# models_dict = Dict{String, Any}()


# for sheet in XLSX.sheetnames(xf)
#     df = DataFrame(XLSX.readtable(data_file, sheet))
#     if sheet == "Bus"
#         bus = Bus{Float64}(length(df.idx))
#         bus.idx = Int32.(df.idx)
#         bus.v = Float64.(df.v0)
#         bus.theta = Float64.(df.a0)
#         bus.vd = zeros(Float64, length(df.idx))
#         bus.vq = zeros(Float64, length(df.idx))
#         bus.i_d = zeros(Float64, length(df.idx))
#         bus.i_q = zeros(Float64, length(df.idx))
#         models_dict["bus"] = bus
#     end
#     if sheet == "Line"
#         line = Line{Float64}(length(df.idx))
#         line.idx = Int32.(df.idx)
#         line.bus1_idx = Int32.(df.bus1)
#         line.bus2_idx = Int32.(df.bus2)
#         line.R = Float64.(df.r)
#         line.X = Float64.(df.x)
#         line.i_d = zeros(Float64, length(df.idx))
#         line.i_q = zeros(Float64, length(df.idx))
#         models_dict["line"] = line
#     end
#     if sheet == "GENCLS"
#         generator = Generator{Float64}(length(df.idx))
#         generator.idx = Int32.(df.idx)
#         generator.bus = Int32.(df.bus)
#         generator.delta = ones(Float64, length(df.idx))
#         generator.omega = ones(Float64, length(df.idx))
#         generator.M = Float64.(df.M)
#         # generator.p_m = Float64.(df.p_m)
#         generator.i_d = zeros(Float64, length(df.idx))
#         generator.i_q = zeros(Float64, length(df.idx))
#         generator.x_d_prime = Float64.(df.xd1)
#         generator.e_q_prime = ones(Float64, length(df.idx))

#         df_pv = DataFrame(XLSX.readtable(data_file, "PV"))
#         generator.p_m = Float64.(df_pv.p0)

#         models_dict["generator"] = generator
#     end

#     if sheet == "Fault"
#         fault = Fault{Float64}(length(df.idx))
#         fault.bus = Int32.(df.bus)
#         fault.r_s = Float64.(df.rf)
#         fault.l_s = Float64.(df.xf)

#         models_dict["fault"] = fault
#     end
#     if sheet == "PQ"
#         load = Load{Float64}(length(df.idx))
#         load.bus = Int32.(df.bus)
#         load.p = Float64.(df.p0)
#         load.q = Float64.(df.q0)

#         models_dict["load"] = load
#     end

#     if sheet == "Slack"
#         slack = Slack(length(df.idx))
#         slack.bus = Int32.(df.bus)

#         models_dict["slack"] = slack
#     end
# end

end #module