#=
    Run TPS simulations
=#

include("KMC.jl")

function propose(model, traj, configObservers, portion, rate=1.0)
    # How to cut the trajectory
    splitTime = rand(Float64)*traj.maxTime
    portion = rand(Bool) # 0 first 1 second
    #portion = false
    splitTime = portion ? splitTime*rate : traj.maxTime-rate*splitTime

    # Split the trajectory and just keep one
    pTraj1, pTraj2 = splitTrajectory(traj, splitTime)
    partialTraj1 = portion ? pTraj2 : pTraj1

    # Propose a new partial trajectory
    pMaxTime = portion ? splitTime : traj.maxTime - splitTime
    pInitial = portion ? partialTraj1.states[end] : partialTraj1.initial
    partialTraj2 = simulation(model, pInitial, pMaxTime, configObservers)
    partialTraj2 = portion ? partialTraj2 : reverseTrajectory(partialTraj2)

    # Join the two trajectories
    newTrajectory = portion ? joinTrajectory(partialTraj1, partialTraj2) :
                              joinTrajectory(partialTraj2, partialTraj1)

    return newTrajectory

end

function gidx(lst, val)
    idx = 1
    for item = lst
        if item < val
            idx += 1
        end
    end

    if idx > size(lst)[1]
        idx = 0
    end

    return idx
end

function splitTrajectory(traj, splitTime::Float64)
    # Find the index where to split
    idx = gidx(traj.times, splitTime)

    # First partial trajectory
    endIdx = idx == 0 ? size(traj.times)[1] : max(1, idx - 1)
    partialInitial1 = copy(traj.initial)
    partialTimes1 = copy(traj.times[1:endIdx])
    partialStates1 = copy(traj.states[1:endIdx])
    partialObservers1 = []
    for observer in traj.observers
        push!(partialObservers1, copy(observer[1:endIdx]))
    end
    partialMaxTime1 = splitTime
    partialTrajectory1 = trajectory(partialInitial1, partialStates1,
                                    partialTimes1, partialMaxTime1,
                                    partialObservers1)

    # Second partial trajectory
    startIdx = idx == 0 ? size(traj.times)[1] : max(1, idx - 1)
    partialInitial2 = copy(traj.states[startIdx])
    partialTimes2 = copy(traj.times[startIdx:end]) .- splitTime
    partialTimes2[1] = 0
    partialStates2 = copy(traj.states[startIdx:end])
    partialObservers2 = []
    for observer in traj.observers
        push!(partialObservers2, copy(observer[startIdx:end]))
    end
    partialMaxTime2 = copy(traj.maxTime) - splitTime
    partialTrajectory2 = trajectory(partialInitial2, partialStates2,
                                    partialTimes2, partialMaxTime2,
                                    partialObservers2)

    return partialTrajectory1, partialTrajectory2
end

function reverseTrajectory(traj)
    states = reverse(traj.states)
    initial = states[1]
    observers = []
    for observer in traj.observers
        push!(observers, reverse(observer))
    end

    times = [0.0]
    sz = size(traj.times)[1]
    for idx in 1:sz-1
        push!(times, traj.maxTime - traj.times[sz+1-idx])
    end

    return trajectory(initial, states, times, traj.maxTime, observers)
end

function joinTrajectory(traj1, traj2)
    initial = copy(traj1.initial)
    maxTime = traj1.maxTime + traj2.maxTime
    if size(traj2.states)[1] > 1
        states = vcat(traj1.states, traj2.states[2:end])
        times = vcat(traj1.times, traj2.times[2:end] .+ traj1.maxTime)
        observers = []
        for i = 1:size(traj1.observers)[1]
            push!(observers, vcat(traj1.observers[i], traj2.observers[i][2:end]))
        end
    else
        states = copy(traj1.states)
        times = copy(traj1.times)
        observers = []
        for i = 1:size(traj1.observers)[1]
            push!(observers, copy(traj1.observers[i]))
        end
    end

    return trajectory(initial, states, times, maxTime, observers)
end


function metropolis(crit1, crit2)
    r = rand(Float64)
    return r < exp(crit2 - crit1) ? true : false
end


function TPS(model, traj, trajMeasures, numSims, maxTime, observers, warmUp=0,
             quiet=false, flipRate = 0.5, ratio = 1)
    # Get initial probability
    crit = criterion(traj, trajMeasures)

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

    acceptance = 0
    portion = true
    for sim = 1:(numSims+warmUp)
        # Propose a new trajectory and calculate observables
        newTraj = propose(model, traj, observers[configObservers], portion, ratio)
        newTrajMeasures = []
        for idx in trajectoryObservers
            observers[idx].currentMeasure = measureObserver(observers[idx], newTraj)
            push!(newTrajMeasures, observers[idx].currentMeasure)
        end
        newCrit = criterion(newTraj, newTrajMeasures)

        # Accept or reject
        accept = metropolis(crit, newCrit)
        if accept
            traj = newTraj
            trajMeasures = newTrajMeasures
            crit = newCrit
            acceptance += 1
        end

        r = rand()
        if r < flipRate
            portion = portion ? false : true
        end

        # Measure observers
        if sim > warmUp
            i = 1
            for idx in configObservers
                measure = timeIntegrate(traj.times, traj.observers[i], maxTime)
                storeMeasure!(observers[idx], measure)
                updateObserver!(observers[idx], measure)
                i += 1
            end
            i = 1
            for idx in trajectoryObservers
                storeMeasure!(observers[idx], trajMeasures[i])
                updateObserver!(observers[idx], trajMeasures[i])
                i += 1
            end
        end

        if !quiet
            println(string("Simulation ", sim, "/", numSims+warmUp, " completed."))
        end
    end
    if !quiet
        println(string("Acceptance rate ", acceptance/(numSims+warmUp), " completed."))
    end

    return traj, trajMeasures, acceptance
end
