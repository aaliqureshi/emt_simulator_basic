module SlackModel

export Slack

mutable struct Slack
    bus::Vector{Int32}

    function Slack(n::Integer)
        new(Vector{Int32}(undef, n))
    end
end

end # module
