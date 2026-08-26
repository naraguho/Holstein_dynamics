using DelimitedFiles  
using LinearAlgebra
using Random
using Distributions



function Half_filling_sorting_Maximize(Qi   :: Vector{Float64}) #Array{Float64, 1} == Vector{Float64}
    if length(Qi) % 2 == 0
        sorted_indices  = sortperm(Qi, rev=true) # without rev = sortperm is ascending order, with rev = sortperm is decending order
        n               = length(Qi)
        half_n          = div(n, 2)
        output          = zeros(n)
        output[sorted_indices[1:half_n]] .= 1
    end
    return output
end
Q1 = [-1, 1.1, -2, 2, 3, 4]
#If Q1 = [-1, 1, -2, 2, 3, 4], function argument should be Vector{Int64}. Julia is strict with this
#println(Half_filling_sorting_Maximize(Q1))

function ConservativeForce_Hamiltonian(Qi   :: Vector{Float64},
                                       ni   :: Vector{Float64},
                                       kappa:: Float64,
                                       k    :: Float64,
                                       g    :: Float64)
  
    Qi_right    = circshift(Qi, 1)
    Qi_left     = circshift(Qi,-1)
    F_lattice   = kappa * ( Qi_right + Qi_left ) + k * Qi
    F_electron  = g * (ni .- 0.5) 
    return -F_lattice + F_electron
end


function kink_counting_ni(ni    :: Vector{Float64})
    boundary_kink   = (ni[1] == ni[end] ? 1 : 0) #Familiar syntax from gnuplot; For partial fitting used it
    ni_left         = circshift(ni, -1)
    kink_number     = sum(ni[1:end-1] .== ni_left[1:end-1]) + boundary_kink #(...) is element wise comparison
    return kink_number
end
#n1 = [1.0,1.0,1,0,1,0,1,0]
#println(kink_counting_ni(n1))

function nbkink(n1)
    kink_00_n1      = (n1[1:end-1] .== 0) .& (n1[2:end] .== 0) #.& is for element wise logical, && is for scalar value boolean operand
    kink_11_n1      = (n1[1:end-1] .== 1) .& (n1[2:end] .== 1)
    kink_signal_n1  = kink_00_n1 .| kink_11_n1
    kink_signal_n1_indices = findall(kink_signal_n1)
    #kink_signal_n1_indices = [1,10,13,30,33]
    sort!(kink_signal_n1_indices)

    nearby = 0
    for i in 1:length(kink_signal_n1_indices)-1
        if kink_signal_n1_indices[i+1] - kink_signal_n1_indices[i] < 5
            nearby += 1
        end
    end
    return nearby
end


function evolution_1DHolsteinDynamics(Qi      :: Vector{Float64},
                                      Pi      :: Vector{Float64},
                                      ni      :: Vector{Float64},
                                      nS      ::Int64,
                                      m       ::Float64,
                                      kappa   ::Float64,
                                      g       ::Float64,
                                      gamma   ::Float64,
                                      T       ::Float64,
                                      dt      ::Float64,
                                      k,
                                      rng)
    kink_number = zeros(nS)
    nbkink_number = zeros(nS)
    kink_sub = zeros(100)
    for step in 1:nS
        #if step % 100 == 0
        #    kink_sub[Int16(step/100)] = kink_counting_ni(ni)                            
        #end
        Force               = ConservativeForce_Hamiltonian(Qi, ni, kappa, k, g) #Vector
        acceleration        = Force / m
        Qi                  += Pi * dt / m + 0.5 * acceleration * dt^2 #a^2 is square, as like C++
        Force_new           = ConservativeForce_Hamiltonian(Qi, ni, kappa, k, g)
        acceleration_new    = Force_new / m
        Pi                  += 1/2 * m * (acceleration + acceleration_new) * dt

        Pi                  += -gamma * Pi *dt
        Pi                  += rand(rng, Normal(0, sqrt(2 * gamma * m * T * dt)), length(Pi))
        #Pi += randn(rng, length(Pi)) * sqrt(2 * gamma * m * T )
        ni                  = Half_filling_sorting_Maximize(Qi)
        kink_number[step]   = kink_counting_ni(ni)
        nbkink_number[step] = nbkink(ni)
        #println(step)

    end
    #print(rand(Normal(0, sqrt(2*gamma*m*T*dt)), length(Pi)))
    return kink_number
end

####################################################################################################
#####                                                                                          #####
#####                              main starts here                                            #####
#####                                                                                          #####
####################################################################################################
function main(T)
    g       = 2.5
    k       = 1.0
    kappa   = 0.3
    nS      = 10000
    m       = 1.0
    gamma   = 0.16
    dt      = 0.1
    T       = T
    nRun    = 200
    dim     = 50
    #Random.seed!(2024)
    #rng = Xoshiro(1)

    P_initial       = zeros(dim)
    n_initial       = zeros(dim)
    n_checkboard    = [i % 2 == 0 ? 1.0 : 0.0 for i in 1:dim]
    #println(n_checkboard)

    global total    = zeros(nS)
    #global Q = zeros(dim)
    for run in 1:nRun
        rng = Xoshiro(2024 + run)
        Q_initial = rand(rng, Normal(0.0, 1.0), dim)
        #print(Q_initial)
        kink_number = evolution_1DHolsteinDynamics(Q_initial, P_initial, n_checkboard, nS, m, kappa, g, gamma, T, dt, k,rng)
        global total += kink_number
        println("($run) is done")
    end
    total /= nRun


    writedlm("san2$T.txt", total)
end

for T in 0.1:0.1:0.5
    main(T)
end