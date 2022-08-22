#=
    The template for 2-level systems
=#

mutable struct Spin2DOriginal{}
    size::Int
    c::Float64
    s::Float64
    state::Array{Bool}
    transitionRates::Array{Float64}
    originalTransitionRates::Array{Float64}
    exp::Float64
end

"""
    Spin(N, c)

Return a two-level spin structure with N sites and c bias.
"""
function Spin2DOriginal(N, c, s)
    return Spin2DOriginal(N, c, s, zeros(N, N), zeros(N^2), zeros(N^2), exp(-s))
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
