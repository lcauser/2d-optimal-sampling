#=
    The template for 2-level systems
=#

mutable struct SEP2DDoob{}
    size::Int
    s::Float64
    state::Array{Bool}
    transitionRates::Array{Float64}
    originalTransitionRates::Array{Float64}
    psi::GPEPS
    sh::Sitetypes
    chi::Int
end

"""
    SEP2D(N)

Return a 2D square lettice for SEP.
"""
function SEP2DDoob(N, s, psi, sh, chi)
    return SEP2DDoob(N, s, zeros(Bool, N, N), zeros(Float64, 2*N*(N+1)), zeros(Float64, 2*N*(N+1)), psi, sh, chi)
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
function transition!(template, idx::Int)
    N = template.size
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

    if x1 > 0 && x1 <= N && y1 > 0 && y1 <= N
        template.state[y1, x1] = 1 - template.state[y1, x1]
    end
    if x2 > 0 && x2 <= N && y2 > 0 && y2 <= N
        template.state[y2, x2] = 1 - template.state[y2, x2]
    end
end

function updateTransitionRates!(template, idx, idx2)
    # Create operator list
    ops = OpList2d(template.size)

    # Find the indexs which can flip
    transitionRates = copy(template.originalTransitionRates)
    for i = 1:length(transitionRates)
        if transitionRates[i] != 0
            # Find the site
            N = template.size
            if i <= N*(N+1)
                y1 = Int(floor((i - 1) / (N+1))) + 1
                y2 = y1
                x1 = i - ((y1-1) * (N+1)) - 1
                x2 = x1 + 1
                dir = false
            elseif i > N*(N+1)
                i2 = i - N*(N+1)
                x1 = Int(floor((i2 - 1) / (N+1))) + 1
                x2 = x1
                y1 = i2 - ((x1-1) * (N+1)) - 1
                y2 = y1 + 1
                dir = true
            end

            # Add the relevant transitions
            if y1 > 0 && y1 <= N && x1 > 0 && x1 <= N && y2 > 0 && y2 <= N && x2 > 0 && x2 <= N
                add!(ops, ["x", "x"], [y1, x1], dir, 1)
                add!(ops, ["id", "id"], [y1, x1], dir, 1)
            elseif y2 > 0 && y2 <= N && x2 > 0 && x2 <= N
                add!(ops, ["x"], [y2, x2], dir, 1)
                add!(ops, ["id"], [y2, x2], dir, 1)
            else
                add!(ops, ["x"], [y1, x1], dir, 1)
                add!(ops, ["id"], [y1, x1], dir, 1)
            end
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
            transitionRates[i] = abs(exp(-template.s) * overlaps[2*j-1] / overlaps[2*j] * transitionRates[i])
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
