#=
    The observer struct is the base for observers
=#

mutable struct observer{}
    name::String
    type::String
    func
    numSims::Int
    measure
    currentMeasure
end

function measureObserver(observer, model)
    return observer.func(model)
end

function measureObserver(observer, traj::trajectory)
    return observer.func(traj)
end

function updateObserver!(observer, measure)
    observer.numSims += 1
    observer.measure += measure
end

function storeMeasure!(observer, measure)
    observer.currentMeasure = measure
end

function getMeasure(observer)
    return observer.measure / observer.numSims
end
