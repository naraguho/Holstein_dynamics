## check the order parameter from low T to high T

using Random
using CUDA
using Distributions
using LinearAlgebra
using DelimitedFiles
using Statistics

function initialize_lattice_kernel(u, uq)

    idx  = (blockIdx().x - 1) * blockDim().x + threadIdx().x # assigining cuda to each linear element in vector 
    ## π order
    sign = mod(idx, 2) - Float64(0.5) #mod(a,2) means 0 or 1 
    @inbounds uq[idx] = sign * Float64(2.0) * u
    #@inbound accelerate programming ; don't check the bound of uq, just follow
    # RHS means just rescailing.

    return nothing
end
# Role of this function : Initialize the checkboard type lattice + rescailing


function initializeHamiltonian_kernel(dim, t0, g, uq, Hmn)

    idx  = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    ## compute the index for left and right neighbors
    lIdx = mod(idx - 2, dim) + 1
    rIdx = mod(idx,     dim) + 1

    @inbounds Hmn[idx, lIdx] = -t0
    @inbounds Hmn[idx, rIdx] = -t0

    @inbounds Hmn[idx, idx] = -g * uq[idx]

    return nothing
end

function updateHamiltonian_kernel(g, uq, Hmn)

    idx  = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    @inbounds Hmn[idx, idx] = -g * uq[idx]
     #ni ~ a_i(dagger, ai), t~aidagger aj

    return nothing
end

function computeForceLattice_kernel(dim, k0, κ, uq, fq) # Onsite force only the lattice
    #Thermalize ; evolving without electron interaction --> reached thermal equilibrium with environment

    idx  = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    ## compute the index for left and right neighbors
    lIdx = mod(idx - 2, dim) + 1
    rIdx = mod(idx,     dim) + 1
    ## compute the force
    @inbounds fq[idx] = -k0 * uq[idx] - κ * (uq[lIdx] + uq[rIdx]) #from the H

    return nothing
end

function computeForce_kernel(dim, fermiF, k0, κ, g, uq, fq, rho, op) #after termalizing

    idx  = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    ## compute the index for left and right neighbors
    lIdx = mod(idx - 2, dim) + 1
    rIdx = mod(idx,     dim) + 1

    nq = Float64(0.0) #computing the electron onsite occupation number
    for i = 1 : dim
        @inbounds nq += rho[idx, i] * rho[idx,  i] * fermiF[i] #fermi factor is 1/(e-betaE + 1)
    end

    @inbounds fq[idx] = -k0 * uq[idx] - κ * (uq[lIdx] + uq[rIdx]) + g * (nq - 0.5) #computing force

    ## compute CDW order parameter
    sign = mod(idx, 2) - Float64(0.5) 
    @inbounds op[idx] = Float64(2.0) * sign * (nq - 0.5)
    #op - order parameter

    return nothing
end

function computeEnergy_kernel(dim, k0, κ, g, uq, vq, ek, ep) #lattice

    idx  = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    ## compute the index for left and right neighbors
    lIdx = mod(idx - 2, dim) + 1

    ## output kinetic energy
    @inbounds ek[idx] = Float64(0.5) * (vq[idx] * vq[idx])
    ## output potential energy
    @inbounds ep[idx] = (uq[idx] * uq[lIdx]) * κ + Float64(0.5) * ((uq[idx] * uq[idx]) * k0 + uq[idx] * g)

    return nothing
end

function computeCDW_kernel(dim, fermiF, rho, nn) #CDW order paramter

    idx  = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    temp = Float64(0.0)
    for i = 1 : dim
        @inbounds temp += rho[idx, i] * rho[idx,  i] * fermiF[i]
    end

    ## assuming zero temperature the electron levels are half filled
    # for i = 1 : Int64(dimSq / 2)
    #     @inbounds temp += rho[idx, i] * rho[idx,  i]
    # end

    nn[idx] = temp

    return nothing
end

function initLatticeUq_kernel(dim, dimSq, sc, ux, uy, uq) #Look up table sine cos, since it's used frequently.

    idx  = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    iIdx = div(idx  - 1, dim) + 1
    jIdx = mod(idx  - 1, dim) + 1

    uxTemp = Float64(0.0)
    uyTemp = Float64(0.0)
    for q = 1 : Int64(dim/2)
        phase = (iIdx + jIdx) * q * sc
        uxTemp += uq[q] * cos(phase + mod(q, 2) * pi * 0.5 )
        uyTemp += uq[q] * cos(phase - mod(q, 2) * pi * 0.5 )
    end
    ux[idx] = uxTemp / dimSq
    uy[idx] = uyTemp / dimSq

    return nothing
end

# function initLattice(dim   :: Int64,
#                      dimSq :: Int64,
#                      eta   :: Float64,
#                      ux    :: CuArray{Float64},
#                      uy    :: CuArray{Float64})

#     uq = Array{Float64}(undef, Int64(dim/2))
#     sc = Float64(2.0 * pi) / dim

#     for kpt = 1 : Int64(dim/2)
#         fa      = kpt / dim
#         uq[kpt] = (eta * eta) / (fa * fa + eta * eta)  +  (eta * eta) / ( (fa- Float64(0.5)) * (fa- Float64(0.5)) + eta * eta)
#     end
#     uNorm = uq' * uq;
#     uq    /= sqrt(uNorm)
#     cu_uq = cu(uq)
#     @cuda threads = nT blocks = nB (
#         initLatticeUq_kernel(dim, dimSq, sc, ux, uy, cu_uq))
#     return nothing
# end



#Evolve sys without electron
function thermalizeSystem(dim   :: Int64,
                          tp    :: Float64,             # temperature of the system
                          gm    :: Float64,             # damping factor for Langevin dynamics
                          k0    :: Float64,             # on-site harmonic potential
                          κ     :: Float64,             # lattice coupling
                          uq    :: CuArray{Float64}, 
                          vq    :: CuArray{Float64},
                          fq    :: CuArray{Float64},
                          dt    :: Float64,
                          nS    :: Int64
                          )
    ## velocity Verlet method is used to update dynamics
    ## compute a(0)
    @cuda threads = nT blocks = nB (
        computeForceLattice_kernel(dim, k0, κ, uq, fq))
    ## compute 0.5 x a(0) x dt
    dvq = Float64(0.5) * dt * fq
    for step = 1 : nS
        ## u(t + dt) = u(t) + v(t + 0.5 x dt) x dt
        uq += (vq + dvq) * dt
        ## compute a(t + dt)
        @cuda threads = nT blocks = nB (
            computeForceLattice_kernel(dim, k0, κ, uq, fq))
        ## v(t + dt) = v(t + 0.5 x dt) + 0.5 x a(t + dt) x dt
        dvqNew = Float64(0.5) * dt * fq
        vq +=   dvq +  dvqNew

        ## generate damping and noise for Langevin dynamics  simple calc can be done in cuda
        etaq = CuArray(rand(rng, Normal(0,1), dim))
        alp  = Float64(exp(-gm * dt))
        sig  = Float64(sqrt((1-alp * alp) * tp))

        ## update velocity using Langevin dynamics
        vq = vq * alp + etaq * sig
    end
    return uq, vq
end

#actual evolution
function evolveSystem(dim    :: Int64,
                      tp     :: Float64,             ## temperature of the system
                      gm     :: Float64,             ## damping factor for Langevin dynamics
                      t0     :: Float64,             ## nearest neighbor hopping
                      g      :: Float64,             ## Holstein coupling
                      k0     :: Float64,             ## on-site harmonic potential
                      κ      :: Float64,             ## lattice coupling
                      uq     :: CuArray{Float64},
                      vq     :: CuArray{Float64},
                      fq     :: CuArray{Float64},
                      dt     :: Float64,
                      nS     :: Int64)

    ## prepare the array for energy output
    ek = CuArray{Float64}(undef, dim)
    ep = CuArray{Float64}(undef, dim)
    eksum = Array{Float64}(undef, 0)
    epsum = Array{Float64}(undef, 0)
    eesum = Array{Float64}(undef, 0)
    op    = Array{Float64}(undef, 0)
    ## prepare output for CDW order
    nn = CuArray{Float64}(undef, dim)
    ## thermalize the lattice system
    # uq, vq = thermalizeSystem(dim, tp, gm, k0, κ, uq, vq, fq, dt, 10000)

    ## prepare the system in π order
    # @cuda threads = nT blocks = nB (
    #     initialize_lattice_kernel(Float64(4.55), uq))

    ## evolve the system with the electron-lattice coupling
    ## initialize the electron Hamiltonian
    Hmn = CUDA.zeros(Float64, (dim, dim))
    # fermiF = vcat(CUDA.ones(Float64, Int64(dimSq/2)), CUDA.zeros(Float64, Int64(dimSq/2)))
    ## compute the initial Hamiltonian
    # @cuda threads = nT blocks = nB (
    #     updateHopping_kernel(dim, dimSq, t0, e_ph, ux, uy, tL, tT, lIdxTbl, tIdxTbl))
    @cuda threads = nT blocks = nB (
        initializeHamiltonian_kernel(dim, t0, g, uq, Hmn))
    ## find eigenvalues and eigenvectors
    vals, vecs = eigen(0.5 * (Hmn + Hmn')) #diagonalizing to compute force
    ## compute Fermi energy and Fermi factor
    CUDA.@allowscalar fermiE = Float64(0.5) * (vals[Int64(dim/2)] + vals[Int64(dim/2)+1])
    CUDA.@allowscalar fermiF = 1.0 ./ ( exp.((vals .-  fermiE ) ./ tp) .+ 1.0)

    ## compute a(0)
    @cuda threads = nT blocks = nB (
        computeForce_kernel(dim, fermiF, k0, κ, g, uq, fq, vecs, nn))

    ## compute 0.5 x a(0) x dt
    dvq = Float64(0.5) * dt * fq

    ## compute the energy of the system
    @cuda threads = nT blocks = nB (
        computeEnergy_kernel(dim, k0, κ, g, uq, vq, ek, ep))
    push!(eksum, reduce(+,ek))
    push!(epsum, reduce(+,ep))
    push!(eesum, vals' * fermiF)
    push!(op, reduce(+,nn) / dim)
    open("evolve.txt", "w") do io
        writedlm(io, Array(uq)')
        writedlm(io, Array(vq)')
        writedlm(io, Array(fq)')
        println("initial energy: $((eksum[end] + epsum[end] + 2.0 * eesum[end]) / dim)")
        #println("initial E: $(reduce(+,ei) / dimSq)")
    end
    open("evolve.txt", "a") do io
    for step = 1 : nS

        uq += (vq + dvq) * dt

        ## update Hamiltonian based on current Q location
        @cuda threads = nT blocks = nB (
            updateHamiltonian_kernel(g, uq, Hmn))

        vals, vecs = eigen(0.5 * (Hmn + Hmn'))
        CUDA.@allowscalar fermiE = Float64(0.5) * (vals[Int64(dim/2)] + vals[Int64(dim/2)+1])
        CUDA.@allowscalar fermiF = 1.0 ./ ( exp.((vals .-  fermiE ) ./ tp) .+ 1.0)
        @cuda threads = nT blocks = nB (
            computeForce_kernel(dim, fermiF, k0, κ, g, uq, fq, vecs, nn))

        dvqNew = Float64(0.5) * dt * fq
        vq +=   dvq +  dvqNew
        dvq =  dvqNew

        etaq = CuArray(rand(rng, Normal(0,1), dim))

        alp  = Float64(exp(-gm *  dt))
        sig  = Float64(sqrt((1 - alp * alp) * tp))

        vq = vq * alp + etaq * sig

        if step % 10 ==0
            writedlm(io, Array(uq)')
            writedlm(io, Array(vq)')
            writedlm(io, Array(fq)')
            @cuda threads = nT blocks = nB (
            computeEnergy_kernel(dim, k0, κ, g, uq, vq, ek, ep))
            push!(eksum, reduce(+,ek))
            push!(epsum, reduce(+,ep))
            push!(eesum, vals' * fermiF)
            push!(op, reduce(+,nn) / dim)
            println("Step: $step done, E: $((eksum[end] + epsum[end] + 2.0 * eesum[end]) / dim)")
        end
    end
    end
    @cuda threads = nT blocks = nB (
            computeCDW_kernel(dim, fermiF, vecs, nn))
    writedlm("CDW.txt", Array(nn)')
    return eksum, epsum, eesum, op
end

###################################################################################################
######
######  main starts here
######
###################################################################################################
# CUDA.allowscalar(true)

dim    = Int64(1024)
t0     = Float64(10.0)
g      = Float64(3.5)
k0     = Float64(1.0)
κ      = Float64(0.18)
dt     = Float64(1e-2)
nS     = Int64(1000)
rng    = Xoshiro(1)

# tp : temperature is in terms of kT
tp = parse(Float64, ARGS[1])
gm = Float64(1.0)

nT = min(512, dim)
nB = cld(dim, nT)

## initialize the system at zero position
uq = CUDA.zeros(Float64, dim)
## initialize the system at randomized position
# uq = CuArray((rand(rng, Float64, dim) .- 0.5) * 1e-5)
## initialize the system from reading a config file
# rawConfig = readdlm("state.txt")
# copyto!(uq, rawConfig[2001,:])

#initLattice(dim, dimSq, 0.1, ux, uy)

## initialize force array
fq = CuArray{Float64}(undef, dim)
## initialize the velocity between -0.5 and 0.5
vq = CuArray(rand(rng, Float64, dim) .- 0.5)
## initialize the velocity according to Maxwell distribution
# vq = CuArray(rand(rng, Rayleigh(sqrt(tp)), dimSq))

## evolve the system
@time ek, ep, ee, op = evolveSystem(dim, tp, gm, t0, g, k0, κ, uq, vq, fq, dt, nS)

open("energy.txt", "w") do io
    writedlm(io, [ek ep ee op])
end
