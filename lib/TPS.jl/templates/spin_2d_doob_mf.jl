#=
    The template for 2-level systems
=#

mutable struct Spin2DDoobMF{}
    size::Int
    c::Float64
    s::Float64
    state::Array{Bool}
    transitionRates::Array{Float64}
    originalTransitionRates::Array{Float64}
    psi::GPEPS
    leftComponent::Float64
end

"""
    Spin(N, c)

Return a two-level spin structure with N sites and c bias.
"""
function Spin2DDoobMF(N, c, s, psi)
    return Spin2DDoobMF(N, c, s, zeros(N, N), zeros(N^2), zeros(N^2), psi, 0)
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
    epsilon = 1e-2
    # Create operator list
    ops = OpList2d(template.size)

    # Find the left component
    leftComponents = [1.0]
    leftComponent = [1.0]
    for i = 2:length(template.originalTransitionRates)
        # Find the site
        y = Int(floor((i-1) / template.size) + 1)
        x = Int((i - 1) % template.size + 1)

        # Multiply
        leftZero = abs(template.psi[y, x][1, 1, 1, 1, 2] * sqrt(1-template.c))
        leftOne = abs(template.psi[y, x][1, 1, 1, 1, 1] * sqrt(template.c))
        ratio = leftOne / leftZero
        ratio = ratio < 1e-5 ? 1e-5 : ratio
        leftComponent *= template.state[y, x] ? ratio : 1
        push!(leftComponents, ratio)
    end
    leftComponent = leftComponent[1]

    # Find the indexs which can flip
    transitionRates = copy(template.originalTransitionRates)
    for i = 1:length(transitionRates)
        if transitionRates[i] != 0
            # Find the site
            y = Int(floor((i-1) / template.size) + 1)
            x = Int((i - 1) % template.size + 1)

            # Find the components
            leftZero = abs(template.psi[y, x][1, 1, 1, 1, 2])
            leftOne = abs(template.psi[y, x][1, 1, 1, 1, 1])
            if template.state[y, x]
                leftComponent2 = leftComponent / leftComponents[i]
            else
                leftComponent2 = leftComponent * leftComponents[i]
            end
            ratio = (leftComponent2 + epsilon) / (leftComponent + epsilon)

            # Update the rate
            transitionRates[i] = abs(exp(-template.s) * ratio * transitionRates[i])
        end


    end

    # Update the rates
    template.transitionRates = transitionRates
    template.leftComponent = leftComponent + epsilon
end

function updateTransitionRates!(template, idx)
    updateOriginalTransitionRates!(template, idx)
    updateTransitionRates!(template, 0, 0)
end

function updateTransitionRates!(template)
    updateOriginalTransitionRates!(template)
    updateTransitionRates!(template, 0, 0)
end
