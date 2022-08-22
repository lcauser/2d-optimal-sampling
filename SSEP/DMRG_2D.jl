include("../lib/TensorNetworks.jl/TensorNetworks.jl")
using HDF5

# Model parameters
N = 6 # System size
c = 0.5 # Temperature
alpha = 1.0 # Boundary
s = 1.0 # Biasing value

# Create an efficient site numbering
labels = zeros(Int, N, N)
i = 1
j = 1
dir = 0
for site = 1:N^2
    labels[i, j] = site
    if dir == 0 && i == N && j == 1
        global j += 1
        global dir = 1 - dir
    elseif dir == 0 && j == 1
        global i += 1
        global dir = 1 - dir
    elseif dir == 0 && i == N
        global j += 1
        global dir = 1 - dir
    elseif dir == 0
        global i += 1
        global j -= 1
    elseif dir == 1 && i == 1 && j == N
        global i += 1
        global dir = 1 - dir
    elseif dir == 1 && i == 1
        global j += 1
        global dir = 1 - dir
    elseif dir == 1 && j == N
        global i += 1
        global dir = 1 - dir
    else
        global i -= 1
        global j += 1
    end
end

# Create hamiltonian
sh = spinhalf()
H = OpList(N^2)
for i = 1:N
    for j = 1:N
        if i == 1 || i == N
            add!(H, ["s+"], [labels[i, j]], -alpha*sqrt(c*(1-c))*exp(-s))
            add!(H, ["s-"], [labels[i, j]], -alpha*sqrt(c*(1-c))*exp(-s))
            add!(H, ["pu"], [labels[i, j]], alpha*(1-c))
            add!(H, ["pd"], [labels[i, j]], alpha*c)
        end
        if j == 1 || j == N
            add!(H, ["s+"], [labels[i, j]], -alpha*sqrt(c*(1-c))*exp(-s))
            add!(H, ["s-"], [labels[i, j]], -alpha*sqrt(c*(1-c))*exp(-s))
            add!(H, ["pu"], [labels[i, j]], alpha*(1-c))
            add!(H, ["pd"], [labels[i, j]], alpha*c)
        end

        if j < N
            add!(H, ["s-", "s+"], [labels[i, j], labels[i, j+1]], -exp(-s))
            add!(H, ["s+", "s-"], [labels[i, j], labels[i, j+1]], -exp(-s))
            add!(H, ["pu", "pd"], [labels[i, j], labels[i, j+1]], 1)
            add!(H, ["pd", "pu"], [labels[i, j], labels[i, j+1]], 1)
        end

        if i < N
            add!(H, ["s-", "s+"], [labels[i, j], labels[i+1, j]], -exp(-s))
            add!(H, ["s+", "s-"], [labels[i, j], labels[i+1, j]], -exp(-s))
            add!(H, ["pu", "pd"], [labels[i, j], labels[i+1, j]], 1)
            add!(H, ["pd", "pu"], [labels[i, j], labels[i+1, j]], 1)
        end
    end
end
H = MPO(sh, H)

# Create activity
K = OpList(N^2)
for i = 1:N
    for j = 1:N
        if i == 1 || i == N
            add!(K, ["s+"], [labels[i, j]], alpha*sqrt(c*(1-c))*exp(-s))
            add!(K, ["s-"], [labels[i, j]], alpha*sqrt(c*(1-c))*exp(-s))
        end
        if j == 1 || j == N
            add!(K, ["s+"], [labels[i, j]], -alpha*sqrt(c*(1-c))*exp(-s))
            add!(K, ["s-"], [labels[i, j]], -alpha*sqrt(c*(1-c))*exp(-s))
        end

        if j < N
            add!(K, ["s-", "s+"], [labels[i, j], labels[i, j+1]], exp(-s))
            add!(K, ["s+", "s-"], [labels[i, j], labels[i, j+1]], exp(-s))
        end

        if i < N
            add!(K, ["s-", "s+"], [labels[i, j], labels[i+1, j]], exp(-s))
            add!(K, ["s+", "s-"], [labels[i, j], labels[i+1, j]], exp(-s))
        end
    end
end
K = MPO(sh, K)

# Create initial guess
psi = randomMPS(2, N^2, 1)
movecenter!(psi, 1)

# List to measure occupations
oplist = OpList(N^2)
for i = 1:N^2
    add!(oplist, ["pu"], [i])
end

# Do DMRG; interactions are "long ranged" so not cutoff, work up to large D
energylast = inner(psi, H, psi)
occslast = inner(sh, psi, oplist, psi)
for D = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
    # DMRG
    global psi, energy = dmrg(psi, H; maxsweeps=1, cutoff=0, maxdim=D, verbose=1)
    global psi, energy = dmrg(psi, H; maxsweeps=1000, cutoff=0, maxdim=D, minsweeps=1, nsites=1, tol=1e-7)

    # Measure occs
    global occs = inner(sh, psi, oplist, psi)

    # Convergence
    converge = 2*abs((energy - energylast) / (energy + energylast)) < 1e-7 ? true : false
    converge = maximum(2*abs.(occs - occslast)) < 1e-6 ? converge : false
    converge && break
    global occslast = occs
    global energylast = energy
end

# Measure observables
activity = inner(psi, K, psi)
expectations = inner(sh, psi, oplist, psi)
occupations = expectations[1:N^2]
occs = zeros(Float64, (N, N))
for i = 1:N
    for j = 1:N
        occs[i, j] = occupations[labels[i, j]]
    end
end