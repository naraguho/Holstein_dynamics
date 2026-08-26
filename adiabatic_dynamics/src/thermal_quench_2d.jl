using LinearAlgebra
using Random
using DelimitedFiles

################################################################################
# Dimensionless 2D semiclassical Holstein model with Langevin bath
#
# Dimensionless lattice equation:
#
#   dQ_i/dt = P_i
#
#   dP_i/dt = (n_i - 1/2) - Q_i - gamma_tilde * P_i + noise
#
# Here Q_i, P_i, and t are dimensionless:
#
#   Q_i = Q_i^phys / Q0,       Q0 = g / K
#   P_i = P_i^phys / P0,       P0 = m Omega Q0
#   t   = Omega * t_phys,      Omega = sqrt(K/m)
#
# Temperature convention:
#
#   T_tilde = k_B T / epsilon
#
# where
#
#   epsilon = g^2 / K
#
# Therefore the Langevin noise satisfies
#
#   <xi_i(t) xi_j(t')>
#   =
#   2 gamma_tilde T_tilde delta_ij delta(t - t')
#
# in dimensionless time.
#
# Electronic Hamiltonian is written in units of t_nn:
#
#   H/t_nn =
#       - sum_<ij> (c_i^† c_j + h.c.)
#       - ep_coeff * sum_i (n_i - 1/2) Q_i
#
# Convention:
#
#   lambda = g^2 / (K W) = epsilon / W
#
# Here we choose
#
#   W = 4 t_nn
#
# so
#
#   epsilon / t_nn = 4 lambda
#
# and therefore
#
#   ep_coeff = g^2 / (K t_nn) = epsilon / t_nn = 4 lambda
#
# The electronic eigenvalues are in units of t_nn, so the electronic Fermi
# temperature must be converted as
#
#   T_el / t_nn = 4 lambda * T_tilde
#
# Important:
#
#   T_tilde is NOT T/t_nn.
#   T_tilde is T/epsilon.
#
################################################################################

#------------------------------------------------------------------------------
# Indexing utilities
#------------------------------------------------------------------------------

@inline function site_index(x::Int, y::Int, Lx::Int, Ly::Int)
    return x + (y - 1) * Lx
end

function build_neighbors_2d(Lx::Int, Ly::Int)
    N = Lx * Ly

    lNbr = zeros(Int, N)
    rNbr = zeros(Int, N)
    uNbr = zeros(Int, N)
    dNbr = zeros(Int, N)

    for y in 1:Ly
        yp = (y == Ly) ? 1 : y + 1
        ym = (y == 1 ) ? Ly : y - 1

        for x in 1:Lx
            xp = (x == Lx) ? 1 : x + 1
            xm = (x == 1 ) ? Lx : x - 1

            i  = site_index(x,  y,  Lx, Ly)
            li = site_index(xm, y,  Lx, Ly)
            ri = site_index(xp, y,  Lx, Ly)
            ui = site_index(x,  yp, Lx, Ly)
            di = site_index(x,  ym, Lx, Ly)

            lNbr[i] = li
            rNbr[i] = ri
            uNbr[i] = ui
            dNbr[i] = di
        end
    end

    return lNbr, rNbr, uNbr, dNbr
end

#------------------------------------------------------------------------------
# Electronic Hamiltonian in units of t_nn
#------------------------------------------------------------------------------

function initializeHamiltonian2D_dimless!(
    N        :: Int,
    ep_coeff :: Float64,
    Qi       :: Vector{Float64},
    Hmn      :: Matrix{Float64},
    lNbr     :: Vector{Int},
    rNbr     :: Vector{Int},
    uNbr     :: Vector{Int},
    dNbr     :: Vector{Int}
)
    fill!(Hmn, 0.0)

    @inbounds for i in 1:N
        Hmn[i, i] = -ep_coeff * Qi[i]

        Hmn[i, lNbr[i]] = -1.0
        Hmn[i, rNbr[i]] = -1.0
        Hmn[i, uNbr[i]] = -1.0
        Hmn[i, dNbr[i]] = -1.0
    end

    return
end

function updateHamiltonian2D_dimless!(
    N        :: Int,
    ep_coeff :: Float64,
    Qi       :: Vector{Float64},
    Hmn      :: Matrix{Float64}
)
    @inbounds for i in 1:N
        Hmn[i, i] = -ep_coeff * Qi[i]
    end

    return
end

#------------------------------------------------------------------------------
# Chemical potential and Fermi function
#------------------------------------------------------------------------------

function total_density_from_mu(vals::Vector{Float64}, mu::Float64, T_el::Float64)
    n = 0.0

    @inbounds for e in vals
        x = (e - mu) / T_el

        if x > 40
            n += 0.0
        elseif x < -40
            n += 1.0
        else
            n += 1.0 / (exp(x) + 1.0)
        end
    end

    return n
end

function find_mu_bisection(
    vals             :: Vector{Float64},
    T_el             :: Float64,
    target_particles :: Float64;
    tol              :: Float64 = 1e-12,
    maxiter          :: Int = 2000
)
    emin = minimum(vals)
    emax = maximum(vals)

    lo = emin - 50.0 * max(T_el, 1.0)
    hi = emax + 50.0 * max(T_el, 1.0)

    f_lo = total_density_from_mu(vals, lo, T_el) - target_particles
    f_hi = total_density_from_mu(vals, hi, T_el) - target_particles

    if f_lo > 0 || f_hi < 0
        error("Bisection bracket failed for chemical potential.")
    end

    for iter in 1:maxiter
        mid = 0.5 * (lo + hi)
        f_mid = total_density_from_mu(vals, mid, T_el) - target_particles

        if abs(f_mid) < tol || abs(hi - lo) < tol
            return mid
        end

        if f_mid > 0
            hi = mid
        else
            lo = mid
        end
    end

    return 0.5 * (lo + hi)
end

function occupations_zeroT(vals::Vector{Float64}, target_particles::Float64)
    N = length(vals)
    occ = zeros(Float64, N)

    n_full = floor(Int, target_particles)
    frac   = target_particles - n_full

    if n_full > 0
        occ[1:n_full] .= 1.0
    end

    if n_full < N && frac > 0
        occ[n_full + 1] = frac
    end

    return occ
end

function fermi_from_mu(vals::Vector{Float64}, mu::Float64, T_el::Float64)
    ff = similar(vals)

    @inbounds for i in eachindex(vals)
        x = (vals[i] - mu) / T_el

        if x > 40
            ff[i] = 0.0
        elseif x < -40
            ff[i] = 1.0
        else
            ff[i] = 1.0 / (exp(x) + 1.0)
        end
    end

    return ff
end

#------------------------------------------------------------------------------
# Occupation and force
#------------------------------------------------------------------------------

function computeOccupation!(
    N      :: Int,
    vecs   :: Matrix{Float64},
    fermiF :: Vector{Float64},
    nn     :: Vector{Float64}
)
    fill!(nn, 0.0)

    @inbounds for m in 1:N
        fm = fermiF[m]

        for i in 1:N
            nn[i] += abs2(vecs[i, m]) * fm
        end
    end

    return
end

function computeForce2D_dimless!(
    N      :: Int,
    vecs   :: Matrix{Float64},
    fermiF :: Vector{Float64},
    Qi     :: Vector{Float64},
    Fi     :: Vector{Float64},
    nn     :: Vector{Float64}
)
    computeOccupation!(N, vecs, fermiF, nn)

    @inbounds for i in 1:N
        Fi[i] = (nn[i] - 0.5) - Qi[i]
    end

    return
end

#------------------------------------------------------------------------------
# Snapshot saving
#------------------------------------------------------------------------------

function save_snapshot(
    run_dir :: String,
    step    :: Int,
    Qi      :: Vector{Float64},
    Pi      :: Vector{Float64},
    nn      :: Vector{Float64},
    Lx      :: Int,
    Ly      :: Int
)
    qfile = joinpath(run_dir, "Q_step_$(lpad(step, 7, '0')).txt")
    pfile = joinpath(run_dir, "P_step_$(lpad(step, 7, '0')).txt")
    nfile = joinpath(run_dir, "n_step_$(lpad(step, 7, '0')).txt")

    writedlm(qfile, reshape(Qi, Lx, Ly))
    writedlm(pfile, reshape(Pi, Lx, Ly))
    writedlm(nfile, reshape(nn, Lx, Ly))

    return
end

#------------------------------------------------------------------------------
# Optional one electronic solve at fixed Q
#------------------------------------------------------------------------------

function compute_density_from_Q!(
    N                :: Int,
    Qi               :: Vector{Float64},
    Hmn              :: Matrix{Float64},
    nn               :: Vector{Float64},
    ep_coeff         :: Float64,
    T_el_dimless     :: Float64,
    target_particles :: Float64
)
    updateHamiltonian2D_dimless!(N, ep_coeff, Qi, Hmn)

    vals, vecs = eigen(Hermitian(Hmn))

    fermiF =
        if T_el_dimless < 1e-12
            occupations_zeroT(vals, target_particles)
        else
            mu = find_mu_bisection(vals, T_el_dimless, target_particles)
            fermi_from_mu(vals, mu, T_el_dimless)
        end

    computeOccupation!(N, vecs, fermiF, nn)

    return vals, vecs, fermiF
end

#------------------------------------------------------------------------------
# Dimensionless time evolution
#------------------------------------------------------------------------------

function evolveSystem2D_dimless!(
    rng,
    Lx                 :: Int,
    Ly                 :: Int,
    nS                 :: Int,
    dt_tilde           :: Float64,
    gamma_tilde        :: Float64,
    T_tilde            :: Float64,
    lambda             :: Float64,
    Qi                 :: Vector{Float64},
    Pi                 :: Vector{Float64},
    Fi                 :: Vector{Float64},
    lNbr               :: Vector{Int},
    rNbr               :: Vector{Int},
    uNbr               :: Vector{Int},
    dNbr               :: Vector{Int};
    filling            :: Float64 = 0.5,
    save_dir           :: Union{Nothing,String} = nothing,
    save_every         :: Int = 100
)
    N = Lx * Ly

    bandwidth_factor = 4.0
    ep_coeff = bandwidth_factor * lambda
    T_el_dimless = bandwidth_factor * lambda * T_tilde

    target_particles = filling * N

    Hmn = zeros(Float64, N, N)
    nn  = zeros(Float64, N)

    initializeHamiltonian2D_dimless!(
        N,
        ep_coeff,
        Qi,
        Hmn,
        lNbr,
        rNbr,
        uNbr,
        dNbr
    )

    vals, vecs = eigen(Hermitian(Hmn))

    fermiF =
        if T_el_dimless < 1e-12
            occupations_zeroT(vals, target_particles)
        else
            mu = find_mu_bisection(vals, T_el_dimless, target_particles)
            fermi_from_mu(vals, mu, T_el_dimless)
        end

    computeForce2D_dimless!(N, vecs, fermiF, Qi, Fi, nn)

    if save_dir !== nothing
        mkpath(save_dir)
        save_snapshot(save_dir, 0, Qi, Pi, nn, Lx, Ly)
    end

    for step in 1:nS
        # Conservative velocity-Verlet step
        dPi = 0.5 .* dt_tilde .* Fi

        @. Qi = Qi + (Pi + dPi) * dt_tilde

        updateHamiltonian2D_dimless!(N, ep_coeff, Qi, Hmn)

        vals, vecs = eigen(Hermitian(Hmn))

        fermiF =
            if T_el_dimless < 1e-12
                occupations_zeroT(vals, target_particles)
            else
                mu = find_mu_bisection(vals, T_el_dimless, target_particles)
                fermi_from_mu(vals, mu, T_el_dimless)
            end

        computeForce2D_dimless!(N, vecs, fermiF, Qi, Fi, nn)

        dPiNew = 0.5 .* dt_tilde .* Fi

        @. Pi = Pi + dPi + dPiNew

        # Exact scalar Ornstein-Uhlenbeck thermostat
        eta = randn(rng, Float64, N)

        alpha = exp(-gamma_tilde * dt_tilde)
        sigma = sqrt((1.0 - alpha * alpha) * T_tilde)

        @. Pi = Pi * alpha + eta * sigma

        if save_dir !== nothing && (step % save_every == 0)
            save_snapshot(save_dir, step, Qi, Pi, nn, Lx, Ly)
        end

        if step % 50 == 0
            println("step $(step) done")
        end
    end

    return
end

#------------------------------------------------------------------------------
# Configuration and main
#------------------------------------------------------------------------------

Base.@kwdef struct SimulationConfig
    n_steps           :: Int = 10_000
    n_runs            :: Int = 20
    Lx                :: Int = 30
    Ly                :: Int = 30
    dt_tilde          :: Float64 = 0.01
    filling           :: Float64 = 0.5
    save_every        :: Int = 10
    lambda            :: Float64 = 0.5
    gamma_tilde       :: Float64 = 0.16
    T_tilde           :: Float64 = 0.005
    T_init_tilde      :: Float64 = 1.0
    base_seed         :: Int = 6_081
    output_root       :: String = normpath(joinpath(@__DIR__, "..", "output"))
end

function validate(config::SimulationConfig)
    config.n_steps >= 0 || error("n_steps must be nonnegative")
    config.n_runs > 0 || error("n_runs must be positive")
    config.Lx > 1 || error("Lx must be greater than 1")
    config.Ly > 1 || error("Ly must be greater than 1")
    config.dt_tilde > 0 || error("dt_tilde must be positive")
    0.0 <= config.filling <= 1.0 || error("filling must lie in [0, 1]")
    config.save_every > 0 || error("save_every must be positive")
    config.lambda >= 0 || error("lambda must be nonnegative")
    config.gamma_tilde >= 0 || error("gamma_tilde must be nonnegative")
    config.T_tilde >= 0 || error("T_tilde must be nonnegative")
    config.T_init_tilde >= 0 || error("T_init_tilde must be nonnegative")
    return
end

function main(config::SimulationConfig = SimulationConfig())
    validate(config)
    N = config.Lx * config.Ly

    base_dir = joinpath(
        config.output_root,
        "Lx$(config.Lx)_Ly$(config.Ly)_lambda$(config.lambda)_" *
        "gamma$(config.gamma_tilde)_T$(config.T_tilde)_" *
        "filling$(config.filling)"
    )

    mkpath(base_dir)

    # Save parameter summary
    open(joinpath(base_dir, "parameters.txt"), "w") do io
        println(io, "n_steps = $(config.n_steps)")
        println(io, "n_runs = $(config.n_runs)")
        println(io, "Lx = $(config.Lx)")
        println(io, "Ly = $(config.Ly)")
        println(io, "N = $(N)")
        println(io, "dt_tilde = $(config.dt_tilde)")
        println(io, "filling = $(config.filling)")
        println(io, "save_every = $(config.save_every)")
        println(io, "lambda = $(config.lambda)")
        println(io, "gamma_tilde = $(config.gamma_tilde)")
        println(io, "T_tilde = $(config.T_tilde)")
        println(io, "T_el_dimless = $(4.0 * config.lambda * config.T_tilde)")
        println(io, "T_init_tilde = $(config.T_init_tilde)")
        println(io, "base_seed = $(config.base_seed)")
        println(io, "Expected kinetic temperature <P_i^2> = $(config.T_tilde)")
        println(io, "Expected kinetic energy per site 0.5<P_i^2> = $(0.5 * config.T_tilde)")
        println(io, "Saved files: Q_step_*.txt, P_step_*.txt, n_step_*.txt")
    end

    # ---------------- Neighbors ----------------
    lNbr, rNbr, uNbr, dNbr = build_neighbors_2d(config.Lx, config.Ly)

    for run in 1:config.n_runs
        println("Starting run $(run)")

        run_dir = joinpath(base_dir, "run_$(run)")
        mkpath(run_dir)

        rng = Xoshiro(config.base_seed + run)

        Qi = sqrt(config.T_init_tilde) .* randn(rng, Float64, N)
        Pi = sqrt(config.T_init_tilde) .* randn(rng, Float64, N)

        Fi = zeros(Float64, N)

        writedlm(joinpath(run_dir, "Q_init.txt"), reshape(Qi, config.Lx, config.Ly))
        writedlm(joinpath(run_dir, "P_init.txt"), reshape(Pi, config.Lx, config.Ly))

        evolveSystem2D_dimless!(
            rng,
            config.Lx,
            config.Ly,
            config.n_steps,
            config.dt_tilde,
            config.gamma_tilde,
            config.T_tilde,
            config.lambda,
            Qi,
            Pi,
            Fi,
            lNbr,
            rNbr,
            uNbr,
            dNbr;
            filling    = config.filling,
            save_dir   = run_dir,
            save_every = config.save_every
        )

        println("Run $(run) finished")
    end

    return
end

#------------------------------------------------------------------------------
# Command-line entry point
#------------------------------------------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__
    default_output = normpath(joinpath(@__DIR__, "..", "output"))
    output_root = length(ARGS) >= 1 ? abspath(ARGS[1]) : default_output
    base_seed = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 6_081

    main(SimulationConfig(output_root = output_root, base_seed = base_seed))
end
