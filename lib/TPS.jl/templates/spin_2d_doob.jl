#=
    The template for 2-level systems
=#

mutable struct Spin2DDoob{}
    size::Int
    c::Float64
    s::Float64
    state::Array{Bool}
    transitionRates::Array{Float64}
    originalTransitionRates::Array{Float64}
    psi::GPEPS
    sh::Sitetypes
    chi::Int
end

"""
    Spin(N, c)

Return a two-level spin structure with N sites and c bias.
"""
function Spin2DDoob(N, c, s, psi, sh, chi)
    return Spin2DDoob(N, c, s, zeros(N, N), zeros(N^2), zeros(N^2), psi, sh, chi)
end


"""
    setState!(template, state)

Change the state in the model
"""
function setState!(template, state)
    template.state = copy(state)
end


"""
    transition!(template, idx)

Transition the state in the model according to the unique identifier idx.
"""
function transition!(template, idx)
    y = Int(floor((idx-1) / template.size) + 1)
    x = Int((idx - 1) % template.size + 1)
    template.state[y, x] = 1 - template.state[y, x]
end

function updateTransitionRates!(template, idx, idx2)
    # Create operator list
    ops = OpList2d(template.size)

    # Find the indexs which can flip
    transitionRates = copy(template.originalTransitionRates)
    for i = 1:length(transitionRates)
        if transitionRates[i] != 0
            # Find the site
            y = Int(floor((i-1) / template.size) + 1)
            x = Int((i - 1) % template.size + 1)
            add!(ops, ["x"], [y, x], false, 1)
            add!(ops, ["id"], [y, x], false, 1)
        end
    end

    # Construct current state as PEPS
    states = [[template.state[i, j] == 1 ? "up" : "dn" for j = 1:template.size] for i = 1:template.size]
    states = productPEPS(template.sh, states)

    # Calculate the expectations
    overlaps = inner(template.sh, template.psi, ops, states; chi=template.chi)

    # Calculate the transition rates
    j = 1
    for i = 1:length(transitionRates)
        if transitionRates[i] != 0
            # Find the site
            y = Int(floor((i-1) / template.size) + 1)
            x = Int((i - 1) % template.size + 1)

            # Update the rate
            factor = template.state[y, x] == 0 ? sqrt(template.c / (1-template.c)) : sqrt((1-template.c) / template.c)
            transitionRates[i] = abs(exp(-template.s) * overlaps[2*j-1] / overlaps[2*j] * transitionRates[i] * factor)
            j += 1
        end
    end

    # Update the rates
    template.transitionRates = transitionRates
end

function updateTransitionRates!(template, idx)
    updateOriginalTransitionRates!(template, idx)
    updateTransitionRates!(template, 0, 0)
end

function updateTransitionRates!(template)
    updateOriginalTransitionRates!(template)
    updateTransitionRates!(template, 0, 0)
end
