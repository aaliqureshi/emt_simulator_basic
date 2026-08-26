module JsonRW

using JSON
using Dates

export save_results, load_results

const SCHEMA_VERSION = 1

# JSON.jl serializes Matrix in column-major order, so a matrix of shape
# (n_var × n_time) would be written as columns (i.e. one inner array per
# time step). To preserve the intended (variable-major) layout on disk,
# we hand JSON a Vector{Vector{Float64}} ourselves.
_rows(m::AbstractMatrix) = [collect(m[i, :]) for i in 1:size(m, 1)]

"""
    save_results(path, sol, sys, address; case_file=nothing, solver=nothing,
                 adaptive=nothing, dt=nothing, extra_metadata=Dict())

Write the simulation result to a pretty-printed JSON file. The state vector is
split by the groups in `address` (`delta`, `omega`, `line_id`, ...) and stored
in `(n_var × n_time)` orientation. A `metadata` block records solver settings
and a `topology` block records bus / line / fault indices so each row of a
state group can be matched back to its physical device.
"""
function save_results(path::AbstractString, sol, sys, address;
                      case_file::Union{Nothing,AbstractString}=nothing,
                      solver=nothing,
                      adaptive::Union{Nothing,Bool}=nothing,
                      dt::Union{Nothing,Real}=nothing,
                      get_states_func=nothing,
                      extra_metadata::AbstractDict=Dict{String,Any}())

    states = get_states_func === nothing ? _states_from_sol(sol) : get_states_func(sol)

    state_groups = Dict{String,Any}()
    for (key, range) in address
        isempty(range) && continue
        state_groups[key] = _rows(states[collect(range), :])
    end

    metadata = Dict{String,Any}(
        "case_file"  => case_file,
        "solver"     => solver === nothing ? nothing : string(solver),
        "adaptive"   => adaptive,
        "dt"         => dt,
        "tspan"      => [first(sol.time), last(sol.time)],
        "n_steps"    => length(sol.time),
        "n_states"   => size(states, 1),
        "t_fault"    => collect(sys.models.fault.t_fault),
        "t_clear"    => collect(sys.models.fault.t_clear),
        "saved_at"   => string(now()),
    )
    for (k, v) in extra_metadata
        metadata[String(k)] = v
    end

    orig = sys.models.bus.orig_idx
    topology = Dict{String,Any}(
        "all_buses"       => collect(Int.(orig)),
        "generator_buses" => collect(Int.(orig[sys.models.generator.bus])),
        "non_slack_buses" => collect(Int.(orig[sys.non_slack_buses])),
        "fault_buses"     => collect(Int.(orig[sys.models.fault.bus])),
        "line_endpoints"  => [[Int(orig[b1]), Int(orig[b2])] for (b1, b2) in
                              zip(sys.models.line.bus1_idx, sys.models.line.bus2_idx)],
    )

    payload = Dict{String,Any}(
        "schema_version" => SCHEMA_VERSION,
        "metadata"       => metadata,
        "topology"       => topology,
        "time"           => collect(Float64, sol.time),
        "states"         => state_groups,
    )

    parent = dirname(path)
    isempty(parent) || isdir(parent) || mkpath(parent)
    open(path, "w") do io
        JSON.print(io, payload, 2)
    end
    return path
end

# Default extractor: tries to call MyDiffEq.get_states if available; otherwise
# expects sol.u to be a Vector of state vectors (one per saved time).
function _states_from_sol(sol)
    if isdefined(Main, :MyDiffEq) && hasproperty(Main.MyDiffEq, :get_states)
        return Main.MyDiffEq.get_states(sol)
    elseif hasproperty(sol, :u)
        u = sol.u
        n_time = length(u)
        n_var  = length(u[1])
        m = Matrix{Float64}(undef, n_var, n_time)
        @inbounds for j in 1:n_time, i in 1:n_var
            m[i, j] = u[j][i]
        end
        return m
    else
        error("Cannot extract states from sol; pass `get_states_func=...`")
    end
end

"""
    load_results(path) -> Dict{Symbol,Any}

Read a JSON results file produced by [`save_results`](@ref). The returned
dictionary has Symbol keys; state-group matrices (e.g. `delta`, `omega`,
`line_id`, ...) are promoted to the top level, alongside `:time`,
`:metadata`, `:topology`, and `:schema_version`. Each state group is a
`Matrix{Float64}` of shape `(n_var × n_time)`.

```julia
d = load_results("results/run_001.json")
d[:delta]                # n_gens × n_time
d[:metadata]["solver"]
d[:topology]["generator_buses"]
```
"""
function load_results(path::AbstractString)
    raw = JSON.parsefile(path)

    out = Dict{Symbol,Any}(
        :schema_version => raw["schema_version"],
        :metadata       => raw["metadata"],
        :topology       => raw["topology"],
        :time           => Vector{Float64}(raw["time"]),
    )

    states = raw["states"]
    for (k, v) in states
        # v is Vector{Vector{Any}} from the JSON parser.
        n_var  = length(v)
        n_time = n_var == 0 ? 0 : length(v[1])
        m = Matrix{Float64}(undef, n_var, n_time)
        @inbounds for i in 1:n_var, j in 1:n_time
            m[i, j] = Float64(v[i][j])
        end
        out[Symbol(k)] = m
    end

    return out
end

end # module
