#=
    Code to create Kinetic Monte Carlo simulations.

    Author: Luke Causer
=#


mutable struct trajectory{}
    initial
    states
    times::Array{Float64}
    maxTime::Float64
    observers::Array{Any}
end

include("observer.jl")

function trajectory(model, initial, idxs, times, maxTime, observers)
    configs, realTimes = reconstruct(model, initial, idxs, times)
    return trajectory(initial, configs, realTimes, maxTime, observers)
end

#function trajectory(model, initial, idxs, times)
#    return trajectory(model, initial, idxs, times, [])
#end


function randomChoice(weights)
    cs = cumsum(weights)
    r = rand(Float64) * cs[end]
    idx = 1
    for val in cs
        if val < r
            idx += 1
        end
    end

    return idx
end


function simulation(model, initial, maxTime, observers)
    # Set the initial state and variables
    setState!(model, initial)
    updateTransitionRates!(model)
    time = 0
    idxs = []
    times = []
    measures = [[] for i=1:size(observers)[1]]

    # Find initial measures
    for i = 1:size(observers)[1]
        push!(measures[i], measureObserver(observers[i], model))
    end

    # Loop through until the max time has passed
    while time < maxTime
        # Find a configuration to transition into and what time
        idx = randomChoice(model.transitionRates)
        escapeRate = sum(model.transitionRates)
        transitionTime = -log(rand(Float64)) / escapeRate
        time += transitionTime

        if time < maxTime
            # Do the transition
            transition!(model, idx)
            updateTransitionRates!(model, idx)

            # Store information
            push!(idxs, idx)
            push!(times, transitionTime)

            # Find measures
            for i = 1:size(observers)[1]
                push!(measures[i], measureObserver(observers[i], model))
            end
        end
    end

    # Store the measures
    for i = 1:size(observers)[1]
        storeMeasure!(observers[i], measures[i])
    end

    return trajectory(model, initial, idxs, times, maxTime, measures)
end

function simulation(model, initial, maxTime)
    return simulation(model, initial, maxTime, [])
end


function reconstruct(model, initial, idxs, times)
    # Preallocate memory
    sz = size(idxs)[1]
    realTimes = Array{Float64}(undef, sz+1)
    states = Array{typeof(model.state)}(undef, sz+1)

    # Set the initial state and variables
    time = 0
    setState!(model, initial)

    # Set starting information
    realTimes[1] = time
    states[1] = initial

    # Loop through changing the configuration
    for i = 1:sz
        transition!(model, idxs[i])
        time += times[i]
        realTimes[i+1] = time
        states[i+1] = copy(model.state)
    end

    return states, realTimes

end


function timeIntegrate(times, measure, maxTime)
    sz = size(measure)[1]
    newTimes = copy(times)
    push!(newTimes, maxTime)
    integral = measure[1]*(newTimes[2]-newTimes[1])
    for i = 2:sz
        integral += measure[i]*(newTimes[i+1]-newTimes[i])
    end

    return integral
end


function KMC(model, numSims, maxTime, observers, quiet=false)
    # Split observers into configuration based and trajectory based
    configObservers = []
    trajectoryObservers = []
    for i = 1:size(observers)[1]
        if observers[i].type == "configuration"
            push!(configObservers, i)
        else
            push!(trajectoryObservers, i)
        end
    end

    # Loop through all simulations
    for sim = 1:numSims
        # Perform a simulation
        initial = initialState(model)
        trajectory = simulation(model, initial, maxTime, observers[configObservers])

        # Measure observables
        i = 1
        for idx in configObservers
            measure = timeIntegrate(trajectory.times, observers[idx].currentMeasure, maxTime)
            updateObserver!(observers[idx], measure)
            i += 1
        end
        for idx in trajectoryObservers
            observers[idx].currentMeasure = measureObserver(observers[idx], trajectory)
            updateObserver!(observers[idx], observers[idx].currentMeasure)
        end

        if !quiet
            println(string("Simulation ", sim, "/", numSims, " completed."))
        end
    end
end
