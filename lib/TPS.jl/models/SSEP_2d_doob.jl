#=
    This creates the type for a 2-level spin system with N sites.
=#

include("../templates/SEP_2d_doob.jl")

"""
    updateOriginalTransitionRates!(template, [idx,])

Updates the transition rates for the set state, only updating local rates
if an identifier idx is provided.
"""
function updateOriginalTransitionRates!(template, idxs::Vector{Int})
    N = template.size
    for idx = idxs
        # Determine x, y coords
        if idx <= N*(N+1)
            y1 = Int(floor((idx - 1) / (N+1))) + 1
            y2 = y1
            x1 = idx - ((y1-1) * (N+1)) - 1
            x2 = x1 + 1
        elseif idx > N*(N+1)
            idx2 = idx - N*(N+1)
            x1 = Int(floor((idx2 - 1) / (N+1))) + 1
            x2 = x1
            y1 = idx2 - ((x1-1) * (N+1)) - 1
            y2 = y1 + 1
        end

        # Check if jump is allowed
        rate::Float64 = 0
        if y1 < 1 || y1 > N || y2 < 1 || y2 > N || x1 < 1 || x1 > N || x2 < 1 || x2 > N
            rate = 0.5
        elseif template.state[y1, x1] == 1 - template.state[y2, x2]
            rate = 1.0
        end

        # Update transition rate
        template.originalTransitionRates[idx] = rate
    end
end

function updateOriginalTransitionRates!(template, idx::Int)
    updateOriginalTransitionRates!(template)
end

function updateOriginalTransitionRates!(template)
    updateOriginalTransitionRates!(template, collect(1:2*template.size*(template.size+1)))
end


"""
    equilibriumState(template)

Generates a configuration from equilibrium.
"""
function equilibriumState(template)
    state = zeros(Bool, template.size, template.size)
    for i = 1:template.size
        for j = 1:template.size
            state[i, j] = rand(Float64) < 0.5 ? 1 : 0
        end
    end
    return state
end
initialState(template) = equilibriumState(template)
