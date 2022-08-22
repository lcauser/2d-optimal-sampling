#=
    Observers for the spin template
=#

# Occupations
function occupationsMeasure(model)
    return convert(Array{Float64}, copy(model.state))
end

function occupationsObserver(N)
    return observer("occupations", "configuration", occupationsMeasure, 0, zeros(N), zeros(N))
end

# Activity
function activityMeasure(traj)
    return size(traj.times)[1]-1
end

function activityObserver()
    return observer("activity", "trajectory", activityMeasure, 0, 0, 0)
end

# Time dependant occupations
function timeOccupationsMeasure(traj, times)
    # Loop through each time
    occupations = []
    for time = times
        # Find the last time smaller or eq to time
        idx = findlast([t <= time for t = traj.times])
        push!(occupations, convert(Array{Float64}, copy(traj.states[idx])))
    end

    return occupations
end

function timeOccupationsObserver(N, times)
    measure(traj) = timeOccupationsMeasure(traj, times)
    initial = []
    for i = 1:size(times)[1]
        push!(initial, zeros(N))
    end
    return observer("time occupations", "trajectory", measure, 0, initial, initial)
end


# Auto correlation function
function autoCorrelationMeasure(traj, times)
    # Loop through each time
    correlations = []
    for time = times
        # Find the last time smaller or eq to time
        idx = findlast([t <= time for t = traj.times])
        push!(correlations, dot(traj.initial, traj.states[idx]))
    end

    return convert(Array{Float64}, correlations)
end

function autoCorrelationObserver(times)
    measure(traj) = autoCorrelationMeasure(traj, times)
    initial = zeros(Float64, size(times))
    return observer("time occupations", "trajectory", measure, 0, initial, initial)
end


# Timescale of AC
function timescaleMeasure(traj)
    timescale = 0.0
    # Loop through each state
    for i = 2:size(traj.states)[1]
        correl = dot(traj.initial, traj.states[i-1])
        timescale += correl * (traj.times[i] - traj.times[i-1])
    end
    correl = dot(traj.initial, traj.states[end])
    timescale += correl * (traj.maxTime - traj.times[end])

    return timescale
end

function timescaleObserver()
    measure(traj) = timescaleMeasure(traj)
    return observer("timescale", "trajectory", measure, 0, 0.0, 0.0)
end

# Measures histograms of activity
function activityHistogramMeasure(traj, times)
    histogram = zeros(size(times)[1])
    for time in traj.times[2:end]
        idx = findfirst([time < t for t = times])
        histogram[idx] += 1
    end
    return histogram
end

function activityHistogramObserver(times)
    measure(traj) = activityHistogramMeasure(traj, times)
    return observer("act_hist", "trajectory", measure, 0, zeros(size(times)[1]),
                    zeros(size(times)[1]))
end
