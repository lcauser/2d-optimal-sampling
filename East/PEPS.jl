using HDF5
include("../lib/TensorNetworks.jl/TensorNetworks.jl")

# System parameters
N = 6 # System size
c = 0.5 # Temperature
s = 1.0 # Biasing
maxdim = 2 # Max PEPS bond dimension
fu_chi = 4 * maxdim^2 # Maximum full update boundary dimension

### Create the Hamiltonian ###
sh = spinhalf()
ops = OpList2d(N)
global kin = 0
for i = 1:N
    for j = 1:N
        p = (i == 1 && j == 1) ? "id" : "n"
        if i <= N-1
            if i == 1 && j == 1
                add!(ops, ["x"], [2, 1], true, sqrt(c*(1-c))*exp(-s))
            else
                add!(ops, [p, "x"], [i, j], true, sqrt(c*(1-c))*exp(-s))
            end
            global kin = kin + 1
        end
        if j <= N-1
            if i == 1 && j == 1
                add!(ops, ["x"], [1, 2], false, sqrt(c*(1-c))*exp(-s))
            else
                add!(ops, [p, "x"], [i, j], false, sqrt(c*(1-c))*exp(-s))
            end
            global kin = kin + 1
        end
    end
end

# Add escape rate terms
global es = 0
for i = 1:N
    for j = 1:N
        p = (i == 1 && j == 1) ? "id" : "n"
        if i <= N-1
            add!(ops, [p, "pu"], [i, j], true, -(1-c))
            add!(ops, [p, "pd"], [i, j], true, -c)
            global es = es + 2
        end
        if j <= N-1
            add!(ops, [p, "pu"], [i, j], false, -(1-c))
            add!(ops, [p, "pd"], [i, j], false, -c)
            global es = es + 2
        end
    end
end

### Simulation ### 
# Create an initial guess
states = [["dn" for i = 1:N] for j = 1:N]
states[1][1] = "up"
states[N][N] = "s"
psi = productPEPS(sh, states)
#psi = randomPEPS(2, N, N, 1)  # Random guess

# Do simple and full updates
psi, energy = simpleupdate(psi, 0.1, sh, ops; maxiter=2000, maxdim=maxdim, saveiter=1000, chi=100, cutoff=1e-8)
psi, energy = simpleupdate(psi, 0.01, sh, ops; maxiter=2000, maxdim=maxdim, saveiter=1000, chi=100, cutoff=1e-8)
psi, energy = fullupdate(psi, ops, 0.01, sh; chi=fu_chi, saveiter=10, maxdim=maxdim, mindim=1, cutoff=1e-6, miniter=1, tol=1e-8, chieval=100)
psi, energy = fullupdate(psi, ops, 0.001, sh; chi=fu_chi, saveiter=10, maxdim=maxdim, mindim=1, cutoff=1e-6, miniter=1, tol=1e-8, chieval=100)


### Measurements ###
# Measure SCGF and activity
observables = inner(sh, psi, ops, psi; maxchi=200)
scgf = sum(observables)
activity = sum(observables[1:kin])

# Measure occupations
ns = OpList2d(N)
for i = 1:N
    for j = 1:N
        add!(ns, ["n"], [i, j], false, 1)
    end
end
ns = inner(sh, psi, ns, psi; maxchi=200)