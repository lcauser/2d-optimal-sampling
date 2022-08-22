using HDF5
include("lib/TensorNetworks.jl/TensorNetworks.jl")
include("lib/TPS.jl/TPS.jl")
include("lib/TPS.jl/models/SSEP_2d_doob.jl")
include("lib/TPS.jl/observers/spin.jl")
include("lib/TPS.jl/observers/doob.jl")

# System parameters
N = 10 # System size
c = 0.5 # temperature
s = 0.01 # Biasing
direct = "D:/SSEP Data/2d/PEPS/alpha = 1.0/c = 0.5/N = 10/s = 1.0.h5" # Where to load the PEPS
num = 10^3 # Number of simulations
chi = 4 # Boundary MPS dimension
maxTime = 50 # Simulation time

# Load in the PEPS
sh = spinhalf()
f = h5open(direct, "r")
psi = read(f, "psi", GPEPS)
activity_peps = read(f, "activity")
close(f)

# Load in the model and observers
model = SEP2DDoob(N, s, psi, sh, chi)
act = activityObserver()
er = escapeRateObserver()
oer = originalEscapeRateObserver()

# Do CTMC simulations
model.state = initialState(model)
acts = []
diffs = []
for i = 1:num
    @time traj = simulation(model, model.state, maxTime, [er, oer])
    act = length(traj.times) - 1
    escapeRate = timeIntegrate(traj.times, traj.observers[1], maxTime)
    originalEscapeRate = timeIntegrate(traj.times, traj.observers[2], maxTime)
    push!(acts, act)
    push!(diffs, originalEscapeRate - escapeRate)
end

# Importance sampling.
diffs = diffs .- minimum(diffs)
k_av = (sum(exp.(-diffs) .*  acts) / sum(exp.(-diffs))) / maxTime