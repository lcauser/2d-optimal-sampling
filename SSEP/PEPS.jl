using HDF5
include("../lib/TensorNetworks.jl/TensorNetworks.jl")

# System parameters
N = 6 # System size
c = 0.5 # Temperature
alpha = 1.0 # Boundary
s = 1.0 # Biasing
maxdim = 2 # Max PEPS bond dimension
fu_chi = 6 * maxdim^2 # Maximum full update boundary dimension

### Create the Hamiltonian ###
sh = spinhalf()
ops = OpList2d(N)
# Add kinetic terms
global kin = 0
for i = 1:N
    for j = 1:N
        if i <= N-1
            add!(ops, ["s+", "s-"], [i, j], true, exp(-s))
            add!(ops, ["s-", "s+"], [i, j], true, exp(-s))
            global kin = kin + 2
        end

        if i == 1 || i == N
            add!(ops, ["s+"], [i, j], true, alpha*sqrt(c*(1-c))*exp(-s))
            add!(ops, ["s-"], [i, j], true, alpha*sqrt(c*(1-c))*exp(-s))
            global kin = kin + 2
        end

        if j <= N-1
            add!(ops, ["s+", "s-"], [i, j], false, exp(-s))
            add!(ops, ["s-", "s+"], [i, j], false, exp(-s))
            global kin = kin + 2
        end
        if j == 1 || j == N
            add!(ops, ["s+"], [i, j], false, alpha*sqrt(c*(1-c))*exp(-s))
            add!(ops, ["s-"], [i, j], false, alpha*sqrt(c*(1-c))*exp(-s))
            global kin = kin + 2
        end
    end
end

# Add escape rate terms
global es = 0
for i = 1:N
    for j = 1:N
        if i <= N-1
            add!(ops, ["pu", "pd"], [i, j], true, -1)
            add!(ops, ["pd", "pu"], [i, j], true, -1)
            global es = es + 2
        end

        if i == 1 || i == N
            add!(ops, ["pu"], [i, j], true, -alpha*(1-c))
            add!(ops, ["pd"], [i, j], true, -alpha*c)
            global es = es + 2
        end

        if j <= N-1
            add!(ops, ["pu", "pd"], [i, j], false, -1)
            add!(ops, ["pd", "pu"], [i, j], false, -1)
            global es = es + 2
        end
        if j == 1 || j == N
            add!(ops, ["pd"], [i, j], false, -alpha*c)
            add!(ops, ["pu"], [i, j], false, -alpha*(1-c))
            global es = es + 2
        end
    end
end

### Simulation ### 
# Create an initial guess
psi = randomPEPS(2, N, N, 1)  # Random guess

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