# data loader
module DataLoader

export load_data

using XLSX, DataFrames, Base.Threads

using ..Models


# normalize a sheet name to a Symbol key
_sym(s::AbstractString) = Symbol(lowercase(s))

function _verify_models(models::NamedTuple) ::Nothing
    """
    Simple test to verify the models for consistency and validity.
    """

    all_buses = models.bus.idx
    gen_buses = models.generator.bus
    load_buses = models.load.bus
    slack_buses = models.slack.bus

    # verify slack bus
    if length(slack_buses) != 1
        error("Oops! There must be exactly one slack bus; got $(length(slack_buses)) !!!")
    end

    ## verify slack bus is not a generator
    if slack_buses[1] in gen_buses
        error("Oops! Slack bus is also a generator - $(slack_buses[1]) !!!")
    end

    ## verify each bus has a classification
    unclassified_buses = setdiff(all_buses, union(gen_buses, load_buses, slack_buses))
    if !isempty(unclassified_buses)
        error("Oops! There are $(length(unclassified_buses)) unclassified buses - $(unclassified_buses) !!!")
    end

    ## ensure only one generator at a bus
    if length(unique(gen_buses)) != length(gen_buses)
        error("Oops! There are more than one generator at the same bus - $(gen_buses) !!!")
    end

    ## ensure only one load at a bus
    if length(unique(load_buses)) != length(load_buses)
        error("Oops! There are more than one load at the same bus - $(load_buses) !!!")
    end

    # TODO: ensure only one line between two same buses

    return nothing
end


function _build_bus(sym_name::Symbol, table::Dict{Symbol, DataFrame})
    """
    Load the bus data from the table.
    """
    df = table[sym_name]
    num_buses = length(df.idx)

    bus = Bus{Float64}(num_buses)
    bus.idx = Int32.(df.idx)
    bus.v = Float64.(df.v0)
    bus.theta = Float64.(df.a0)
    bus.vd = zeros(Float64, length(df.idx))
    bus.vq = zeros(Float64, length(df.idx))

    ## update voltage of PV buses
    df_pv = table[:pv]
    pv_buses = Int32.(df_pv.bus)
    bus.v[pv_buses] = Float64.(df_pv.v0)

    return bus
end

function _build_line(sym_name::Symbol, table::Dict{Symbol, DataFrame})
    """
    Load the line data from the table.
    """
    df = table[sym_name]
    num_lines = length(df.idx)

    line = Line{Float64}(num_lines)
    line.idx = Int32.(collect(1:num_lines))
    line.bus1_idx = Int32.(df.bus1)
    line.bus2_idx = Int32.(df.bus2)
    line.R = max.(Float64.(df.r), 1e-8)
    x_vals = Float64.(df.x)
    line.X = ifelse.(x_vals .== 0.0, 1e-8, x_vals)
    line.B = Float64.(df.b)
    # line.B = max.(Float64.(df.b), 1e-6)
    line.C = line.B ./ (2*pi*60)
    line.L = line.X ./ (2*pi*60)
    line.tap = hasproperty(df, :tap) ? Float64.(df.tap) : ones(Float64, num_lines)
    line.i_d = zeros(Float64, num_lines)
    line.i_q = zeros(Float64, num_lines)

    ## TODO: combine parallel lines

    return line
end

function _build_gencls(sym_name::Symbol, table::Dict{Symbol, DataFrame})
    """
    Load the generator data from the table.
    """
    df = table[sym_name]
    num_gens = length(df.idx)

    generator = Generator{Float64}(num_gens)
    generator.idx = Int32.(collect(1:num_gens))
    generator.bus = Int32.(df.bus)
    generator.delta = ones(Float64, num_gens)
    generator.omega = ones(Float64, num_gens)
    generator.M = Float64.(df.M)
    generator.x_d_prime = Float64.(df.xd1)
    generator.d = Float64.(df.D)
    generator.e_q_prime = ones(Float64, num_gens)
    generator.i_d = zeros(Float64, num_gens)
    generator.i_q = zeros(Float64, num_gens)
    generator.p_m = zeros(Float64, num_gens)
    generator.q_m = zeros(Float64, num_gens)

    ## get p_m and q_m from PV table
    # insert this data as gen_pm and gen_qm
    df_pv = table[:pv]
    pv_buses = Int32.(df_pv.bus)

    pmap = Dict(Int32(b) => Float64(p) for (b,p) in zip(pv_buses, Float64.(df_pv.p0)))
    qmap = Dict(Int32(b) => Float64(q) for (b,q) in zip(pv_buses, Float64.(df_pv.q0)))

    @inbounds for j in 1:num_gens
        b = generator.bus[j]
        generator.p_m[j] = get(pmap, b, 0.0)
        generator.q_m[j] = get(qmap, b, 0.0)
    end

    return generator
end

function _build_fault(sym_name::Symbol, table::Dict{Symbol, DataFrame})
    """
    Load the fault data from the table.
    Note: fault impedance is added a small impedance to avoid division by zero.
    """
    df = table[sym_name]
    num_faults = length(df.idx)

    r_open = [1e10]
    x_open = [1e10]
    Ω = 2*pi*60

    fault = Fault{Float64}(num_faults)
    fault.bus = Int32.(df.bus)

    # fault impedance - add a small impedance to avoid division by zero
    fault.r_fault = max.(Float64.(df.rf), 1e-4)
    fault.x_fault = max.(Float64.(df.xf), 1e-4)
    fault.l_fault = fault.x_fault ./ Ω
    
    fault.t_fault = Float64.(df.tf)
    fault.t_clear = Float64.(df.tc)

    fault.r_open = Float64.(r_open)
    fault.x_open = Float64.(x_open)
    fault.l_open = fault.x_open ./ Ω

    return fault
end

function _build_slack(sym_name::Symbol, table::Dict{Symbol, DataFrame})
    """
    Load the slack data from the table.
    """
    df = table[sym_name]
    num_slack = length(df.idx)

    slack = Slack(num_slack)
    slack.bus = Int32.(df.bus)

    return slack
end

function _build_pq(sym_name::Symbol, table::Dict{Symbol, DataFrame})
    """
    Load the PQ data from the table.
    """
    df_pq = table[sym_name]
    df_bus = table[:bus]
    df_pv = table[:pv]
    df_slack = table[:slack]
    
    all_buses = df_bus.idx
    gen_buses = df_pv.bus
    slack_buses = df_slack.bus

    non_generator_buses = setdiff(all_buses, gen_buses, slack_buses)
    load_buses = Int32.(df_pq.bus)
    connection_buses = setdiff(non_generator_buses, load_buses)

    ## pq_buses are buses with loads and connection buses
    pq_buses = vcat(load_buses, connection_buses)

    ## for connection bus, p and q 
    # loads are zero
    zeros_padded = zeros(Float64, length(connection_buses))

    load = Load{Float64}(length(pq_buses))
    load.bus = Int32.(pq_buses)
    load.p = Float64.(vcat(df_pq.p0, zeros_padded))
    load.q = Float64.(vcat(df_pq.q0, zeros_padded))
    load.y = zeros(Complex{Float64}, length(pq_buses))

    return load
end

function _process_sheet(sym_name::Symbol, table::Dict{Symbol, DataFrame})
    """
    Load Models from the table.
    """
    if sym_name == :bus
        return :bus => _build_bus(sym_name, table)
    elseif sym_name == :line
        return :line => _build_line(sym_name, table)
    elseif sym_name == :gencls
        return :generator => _build_gencls(sym_name, table)
    elseif sym_name == :fault
        return :fault => _build_fault(sym_name, table)
    elseif sym_name == :slack    
        return :slack => _build_slack(sym_name, table)
    elseif sym_name == :pq
        return :load => _build_pq(sym_name, table)
    end
end




# main function
function load_data(xlsx_file_path::String)
    """
    Load data from an Excel file and return a named tuple of models.
    """
    isfile(xlsx_file_path) || error("Case file not found: $xlsx_file_path")

    tables = Dict{Symbol, DataFrame}()
    XLSX.openxlsx(xlsx_file_path) do xf
        sheet_names = XLSX.sheetnames(xf)
        for sheet in sheet_names
            ws = xf[sheet]
            df = DataFrame(XLSX.gettable(ws))
            tables[_sym(sheet)] = df
        end
    end
    
    # TODO: change the work vector as more devices are added!
    ## Only GENCLS name is supported currently!
    work = [:bus, :line, :fault, :slack, :gencls, :pq]
    results = Vector{Pair{Symbol, Any}}(undef, length(work))
    @sync for (i, key) in enumerate(work)
        @spawn begin
            try
                results[i] = _process_sheet(key, tables)
            catch err
                rethrow(ArgumentError("Failed processing sheet '$key': $(err)"))
            end
        end
    end

    models = (; results...)
    #verify
    _verify_models(models)
    return models
end

end #module