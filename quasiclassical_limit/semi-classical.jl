#using Plots: Plot
using LinearAlgebra
using Random
using LinearAlgebra
using DelimitedFiles

#using Plots


########################################################################################
#####                                                                              #####
#####                        Holstein model simulation                             #####
#####                                                                              #####
########################################################################################
function computeForceLattice!(dim, k0, κ, Qi, Fi)
    for idx in 1:dim
        lIdx = mod(idx - 2, dim) + 1
        rIdx = mod(idx,     dim) + 1
        Fi[idx] = -k0 * Qi[idx] - κ * (Qi[lIdx] + Qi[rIdx])
    end
    return
end

function thermalize!(rng,
                     dim   :: Int64,
                     stp   :: Int64,
                     dt    :: Float64,
                     γ     :: Float64,
                     tp    :: Float64,
                     k0    :: Float64,
                     κ     :: Float64,
                     Qi    :: Array{Float64},
                     Pi    :: Array{Float64},
                     Fi    :: Array{Float64})
    computeForceLattice!(dim, k0, κ, Qi, Fi)
    dPi = Float64(0.5) * dt * Fi
    for j in 1:stp
        ## Q(t + dt) = Q(t) + P(t + 0.5 * dt) * dt
        Qi += (Pi + dPi) * dt
        ## compute F(t + dt)
        computeForceLattice!(dim, k0, κ, Qi, Fi)
        ## P(t + dt) = P(t + 0.5 * dt) + 0.5 x F(t + dt) * dt
        dPiNew = Float64(0.5) * dt * Fi
        Pi +=   dPi +  dPiNew
        ## Langevin dynamics
        ## generate damping and noise
        η =  randn(rng, Float64, dim)
        α  = exp(-γ * dt)
        σ  = sqrt((1 - α * α) * tp)
        ## update velocity using Langevin dynamics
        Pi = Pi * α + η * σ
    end
    println("thermalization is run for: $(stp) steps")
    println("current temperature: $((Pi' * Pi) / dim)")
    return
end
########################################################################################
#####                                                                              #####
#####                        Holstein model simulation                             #####
#####                                                                              #####
########################################################################################
function initializeHamiltonian!(dim   :: Int64,
                                t0    :: Float64,
                                g     :: Float64,
                                Qi    :: Array{Float64},
                                Hmn   :: Array{Float64, 2})
    for idx in 1:dim
        lIdx = mod(idx - 2, dim) + 1
        rIdx = mod(idx,     dim) + 1
        Hmn[idx, lIdx] = -t0
        Hmn[idx, rIdx] = -t0
        Hmn[idx, idx]  = -g * Qi[idx]
    end
    return
end

function updateHamiltonian!(dim   :: Int64,
                            g     :: Float64,
                            Qi    :: Array{Float64},
                            Hmn   :: Array{Float64, 2})
    for idx = 1:dim
        Hmn[idx, idx] = -g * Qi[idx]
    end
    return
end

function computeOccupation(dim    :: Int64,
                           idx    :: Int64,
                           ρ      :: Array{Float64, 2},
                           fermiF :: Array{Float64})
    nQ = Float64(0.0)
    for j in 1:dim
        nQ += ρ[idx, j] * ρ[idx,  j] * fermiF[j]
    end
    return nQ
end

function computeForce!(dim    :: Int64,
                       g      :: Float64,
                       k0     :: Float64,
                       κ      :: Float64,
                       ρ      :: Array{Float64, 2},
                       fermiF :: Array{Float64},
                       Qi     :: Array{Float64},
                       Fi     :: Array{Float64})
    ## array for saving electron occupation number
    nn = zeros(Float64, dim)
    for idx in 1:dim
        lIdx = mod(idx - 2, dim) + 1
        rIdx = mod(idx,     dim) + 1
        ## find the occupation number of electrons
        nQ = computeOccupation(dim, idx, ρ, fermiF)
        nn[idx] = nQ
        Fi[idx] = -k0 * Qi[idx] - κ * (Qi[lIdx] + Qi[rIdx]) + g * (nQ - Float64(0.5))
    end
    return nn
end


function countDefect(dim :: Int64,
                      nn :: Array{Float64})
    nDefect = Int64(0)

    #Boundary kink consideration
    if nn[1] > Float64(0.5) && nn[dim] > Float64(0.5)
        nDefect += 1
    elseif nn[1] < Float64(0.5) && nn[dim] < Float64(0.5)
        nDefect += 1
    else
    end
    
    for i in 2:dim
        if nn[i] > Float64(0.5) && nn[i - 1] > Float64(0.5)
            nDefect += 1
        elseif nn[i] < Float64(0.5) && nn[i - 1] < Float64(0.5)
            nDefect += 1
        else
            continue
        end
    end
    return nDefect
end

function evolveSystem!(rng,
                       dim   :: Int64,
                       nS    :: Int64,
                       dt    :: Float64,
                       γ     :: Float64,
                       tp    :: Float64,
                       t0    :: Float64,
                       g     :: Float64,
                       k0    :: Float64,
                       κ     :: Float64,
                       Qi    :: Array{Float64},
                       Pi    :: Array{Float64},
                       Fi    :: Array{Float64})
    ## prepare output array for domain walls
    nD = zeros(Int64, nS + 1)
    ## initialize the Hamiltonian
    Hmn    = zeros(Float64, (dim, dim))
    initializeHamiltonian!(dim, t0, g, Qi, Hmn)
    ## find eigenvalues and eigenvectors
    vals, vecs = eigen(Hermitian(Hmn))
    ## compute Fermi energy and Fermi factor for the half filling
    fermiE = Float64(0.5) * (vals[dim ÷ 2] + vals[dim ÷ 2 + 1])
    fermiF = 1.0 ./ ( exp.((vals .-  fermiE ) ./ tp) .+ 1.0)
    ## compute initial forces when the electrons are present
    nn     = computeForce!(dim, g, k0, κ, vecs, fermiF, Qi, Fi)
    #dPi    = Float64(0.5) * dt * Fi
    nD[1]  = countDefect(dim, nn)
    for j in 1:nS
        ## 0.5 * F(t) * dt
        dPi = Float64(0.5) * dt * Fi
        ## Q(t + dt) = Q(t) + P(t + 0.5 * dt) * dt
        Qi += (Pi + dPi) * dt
        ## compute F(t + dt) from new lattice and electron configuration
        updateHamiltonian!(dim, g, Qi, Hmn)
        vals, vecs = eigen(Hermitian(Hmn))
        fermiE = Float64(0.5) * (vals[dim ÷ 2] + vals[dim ÷ 2 + 1])
        fermiF = 1.0 ./ ( exp.((vals .-  fermiE ) ./ tp) .+ 1.0)
        ## occupation number of electrons is computed
        nn     = computeForce!(dim, g, k0, κ, vecs, fermiF, Qi, Fi) # Fi is updated inside of ftn
        ## P(t + dt) = P(t + 0.5 * dt) + 0.5 x F(t + dt) * dt
        dPiNew = Float64(0.5) * dt * Fi
        Pi    += dPi +  dPiNew
        ## update velocity using Langevin dynamics
        η =  randn(rng, Float64, dim)
        α  = exp(-γ * dt)
        σ  = sqrt((1 - α * α) * tp)
        Pi = Pi * α + η * σ

        nD[j + 1] = countDefect(dim, nn)
        if j % 100 ==0
            println("step $(j) done!")
        end
    end
    return nD
end
########################################################################################
######                                                                           #######
######                             main starts here                              #######
######                                                                           #######
########################################################################################

function main(rng,
              tp  :: Float64)

    ## system info
    ## nS:   evolve sytem for $nS steps
    ## nR:   average over $nR runs
    ## dim:  $dim sites in the 1D system
    nS     = Int64(500000)
    nR     = Int64(10)
    dim    = Int64(1000)
    ## hopping t
    t0     = Float64(0.01)
    ## Holstein coupling
    ## g * (n - 0.5)
    g      = Float64(2.50)
    ## lattice coupling
    ## κ  * Qi * Qj * 0.5
    ## k0 * Qi * Qi * 0.5
    k0     = Float64(1.00)
    κ      = Float64(0.3)
    dt     = Float64(0.02)
    ## Langevin dynamics parameters
    ## damping factor γ
    ## temperature of the system tp
    γ      = Float64(0.16)
    tp     = tp
    ## prepare output averaged
    global tD     = zeros(Float64, nS  + 1)
    for run in 1:nR

        ## two ways to initilize the lattice

        ## (1) randomly initialize the lattice using a normal distribution N(0,1)
        Qi = randn(Float64, dim)
        Pi = randn(Float64, dim)
        Fi = randn(Float64, dim)

        ## (2) thermalize the system to high temperature using Langevin dynamics
        ## set the temperature of the system is high
        # Qi = zeros(Float64, dim)
        # Pi = zeros(Float64, dim)
        # Fi = zeros(Float64, dim)
        # temperature = 1E1
        # thermalize!(rng, dim, 5_000, dt, γ, temperature, k0, κ, Qi, Pi, Fi)

        ## start the simulation for Holstein model
        nD = evolveSystem!(rng, dim, nS, dt, γ, tp, t0, g, k0, κ, Qi, Pi, Fi)
        tD += nD
    end
    tD /= nR
    writedlm("y0.01.txt", tD)
    return
end


slurm_job_id = get(ENV, "SLURM_JOB_ID", "default")
seed = parse(Int, slurm_job_id)
rng = Xoshiro(seed)



main(rng, 0.5)
