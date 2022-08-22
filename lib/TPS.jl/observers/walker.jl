#=
    Observers for the walker template
=#

# position
function positionMeasure(traj, times)
    # Loop through each time
    positions = []
    for time = times
        # Find the last time smaller or eq to time
        idx = findlast([t <= time for t = traj.times])
        push!(positions, convert(Float64, traj.states[idx]))
    end

    return positions
end

function positionObserver(times)
    measure(traj) = positionMeasure(traj, times)
    initial = []
    for i = 1:size(times)[1]
        push!(initial, 0)
    end
    return observer("position", "trajectory", measure, 0, initial, initial)
end


# RMSD
function MSDMeasure(traj, times)
    # Loop through each time
    positions = []
    for time = times
        # Find the last time smaller or eq to time
        idx = findlast([t <= time for t = traj.times])
        push!(positions, convert(Float64, traj.states[idx]^2))
    end

    return positions
end

function MSDObserver(times)
    measure(traj) = MSDMeasure(traj, times)
    initial = []
    for i = 1:size(times)[1]
        push!(initial, 0)
    end
    return observer("MSD", "trajectory", measure, 0, initial, initial)
end
