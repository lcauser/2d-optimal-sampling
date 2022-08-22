#=
    The template for 2-level systems
=#

mutable struct SEP2D{}
    size::Int
    state::Array{Bool}
    transitionRates::Array{Float64}
end

"""
    SEP2D(N)

Return a 2D square lettice for SEP.
"""
function SEP2D(N)
    return SEP2D(N, zeros(Bool, N, N), zeros(Float64, 2*N*(N+1)))
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
