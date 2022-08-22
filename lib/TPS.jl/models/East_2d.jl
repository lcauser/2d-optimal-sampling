#=
    This creates the type for a 2-level spin system with N sites.
=#

include("../templates/spin_2d.jl")

"""
    updateTransitionRates!(template, [idx,])

Updates the transition rates for the set state, only updating local rates
if an identifier idx is provided.
"""
function updateTransitionRates!(template, idxs::Vector{Int})
    for idx = idxs
        # Determine x, y coords
        y = Int(floor((idx-1) / template.size) + 1)
        x = Int((idx - 1) % template.size + 1)

        # Find the kinetic constraint
        constraint::Float64 = 0
        if y > 1
            constraint += template.state[y-1, x]
        end
        if x > 1
            constraint += template.state[y, x-1]
        end

        # Calculate the transition rates
        if template.state[y, x] == 1
            template.transitionRates[idx] = constraint * (1 - template.c)
        else
            template.transitionRates[idx] = constraint * template.c
        end
    end
end

function updateTransitionRates!(template, idx::Int)
    idxs = [idx]

    if idx % template.size != 0
        push!(idxs, idx+1)
    end
    if template.size*(template.size-1) > idx
        push!(idxs, idx + template.size)
    end

    updateTransitionRates!(template, idxs)
end

function updateTransitionRates!(template)
    updateTransitionRates!(template, collect(1:template.size^2))
end


"""
    equilibriumState(template)

Generates a configuration from equilibrium.
"""
function equilibriumState(template)
    state = zeros(Bool, template.size, template.size)
    for i = 1:template.size
        for j = 1:template.size
            state[i, j] = rand(Float64) < template.c ? 1 : 0
        end
    end
    state[1, 1] = 1
    return state
end
initialState(template) = equilibriumState(template)
