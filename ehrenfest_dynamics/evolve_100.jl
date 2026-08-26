using LinearAlgebra
using Random
using DelimitedFiles
using CUDA
using CUDA: memory_status, CuEvent, record, synchronize
using CUDA.CUSPARSE
using SparseArrays


# ============================================================
# Initialization
# ============================================================

function init_current!(dim, Q, P, ρ)
    # Initialize Q and P on GPU
    Q .= 1e-2 .* CUDA.randn(Float64, dim)
    P .= 1e-2 .* CUDA.randn(Float64, dim)

    # Compute ρ on CPU
    tot_k = 2π * (0:dim-1) / dim .- π
    sum_k = filter(k -> -π/2 < k <= π/2, tot_k)

    ρ_host = Array{ComplexF64}(undef, dim, dim)

    for i in 1:dim
        for j in 1:dim
            acc = ComplexF64(0.0, 0.0)
            for k in sum_k
                acc += exp(im * k * (i - j))
            end
            ρ_host[i, j] = acc / dim
        end
    end

    # Transfer ρ to GPU
    copyto!(ρ, ρ_host)
end


function init_gs_deep_quench!(dim, Q, P, ρ; tp=1e-10)
    Q .= 1e-4 .* CUDA.randn(Float64, dim)
    P .= 1e-4 .* CUDA.randn(Float64, dim)

    H0 = hopping_Hamiltonian_1D_dense(dim)

    ρ .= compute_density_matrix_C_cuda(H0, tp, 0.0)
end


function compute_density_matrix_C_cuda(H_gpu, kT, mu)
    H = Array(H_gpu)

    vals, vecs = eigen(H)
    dim = length(vals)

    occ = @. 1 / (exp((vals - mu) / kT) + 1)

    D = zeros(ComplexF64, dim, dim)

    for a in 1:dim
        for b in a:dim
            s = ComplexF64(0.0, 0.0)

            for m in 1:dim
                s += occ[m] * conj(vecs[b, m]) * vecs[a, m]
            end

            D[a, b] = s

            if a != b
                D[b, a] = conj(s)
            end
        end
    end

    return CuArray(D)
end


# ============================================================
# Hamiltonian
# ============================================================

function hopping_Hamiltonian_1D_dense(dim::Int)
    H = zeros(Float64, dim, dim)

    for i in 1:dim
        j_plus  = mod1(i + 1, dim)
        j_minus = mod1(i - 1, dim)

        H[i, j_plus]  = -1.0
        H[i, j_minus] = -1.0
    end

    return CuArray(H)
end


function onsite_Hamiltonian_1D(dim, Q)
    return -Diagonal(Q)
end


# ============================================================
# Equations of motion
# ============================================================

function compute_dQ(Q, P, r, n)
    return r .* P
end


function compute_dP(Q, P, r, n)
    return r .* (n .- 0.5) .- r .* Q
end


# ============================================================
# RK4 time evolution
# ============================================================

function evolve!(
    dim,
    Q,
    P,
    ρ,
    dt,
    r,
    H_hop,
    λ,
    ρ2,
    Kρ_sum,
    Q2,
    KQ_sum,
    P2,
    KP_sum
)

    # ----------------------------
    # RK4 stage 1
    # ----------------------------
    n = real.(diag(ρ))

    KQ = dt .* compute_dQ(Q, P, r, n)
    KP = dt .* compute_dP(Q, P, r, n)

    Kρ = -im * dt .* (
        (H_hop * ρ - ρ * H_hop)
        .+ 4λ .* (reshape(Q, 1, :) .- reshape(Q, :, 1)) .* ρ
    )

    @. ρ2 = ρ + 0.5 * Kρ
    @. Kρ_sum = Kρ / 6

    @. Q2 = Q + 0.5 * KQ
    @. KQ_sum = KQ / 6

    @. P2 = P + 0.5 * KP
    @. KP_sum = KP / 6


    # ----------------------------
    # RK4 stage 2
    # ----------------------------
    n = real.(diag(ρ2))

    KQ = dt .* compute_dQ(Q2, P2, r, n)
    KP = dt .* compute_dP(Q2, P2, r, n)

    Kρ = -im * dt .* (
        (H_hop * ρ2 - ρ2 * H_hop)
        .+ 4λ .* (reshape(Q2, 1, :) .- reshape(Q2, :, 1)) .* ρ2
    )

    @. ρ2 = ρ + 0.5 * Kρ
    @. Kρ_sum += Kρ / 3

    @. Q2 = Q + 0.5 * KQ
    @. KQ_sum += KQ / 3

    @. P2 = P + 0.5 * KP
    @. KP_sum += KP / 3


    # ----------------------------
    # RK4 stage 3
    # ----------------------------
    n = real.(diag(ρ2))

    KQ = dt .* compute_dQ(Q2, P2, r, n)
    KP = dt .* compute_dP(Q2, P2, r, n)

    Kρ = -im * dt .* (
        (H_hop * ρ2 - ρ2 * H_hop)
        .+ 4λ .* (reshape(Q2, 1, :) .- reshape(Q2, :, 1)) .* ρ2
    )

    @. ρ2 = ρ + Kρ
    @. Kρ_sum += Kρ / 3

    @. Q2 = Q + KQ
    @. KQ_sum += KQ / 3

    @. P2 = P + KP
    @. KP_sum += KP / 3


    # ----------------------------
    # RK4 stage 4
    # ----------------------------
    n = real.(diag(ρ2))

    KQ = dt .* compute_dQ(Q2, P2, r, n)
    KP = dt .* compute_dP(Q2, P2, r, n)

    Kρ = -im * dt .* (
        (H_hop * ρ2 - ρ2 * H_hop)
        .+ 4λ .* (reshape(Q2, 1, :) .- reshape(Q2, :, 1)) .* ρ2
    )

    @. Kρ_sum += Kρ / 6
    @. KQ_sum += KQ / 6
    @. KP_sum += KP / 6


    # ----------------------------
    # Update variables
    # ----------------------------
    @. ρ += Kρ_sum
    @. Q += KQ_sum
    @. P += KP_sum
end


# ============================================================
# Simple saving function
# ============================================================

function save_vector_line(filename::String, vec)
    open(filename, "a") do io
        println(io, join(vec, " "))
    end
end


# ============================================================
# Main simulation
# ============================================================

function main(seed::Int, save_base::String)
    Random.seed!(seed)

    dim = 1024
    dt = 0.01
    nS = 1200001
    save_interval = 200

    λ = 0.6
    r = 0.3

    Q = CUDA.zeros(Float64, dim)
    P = CUDA.zeros(Float64, dim)
    ρ = CUDA.zeros(ComplexF64, dim, dim)

    ρ2     = similar(ρ)
    Kρ_sum = similar(ρ)

    Q2     = similar(Q)
    KQ_sum = similar(Q)

    P2     = similar(P)
    KP_sum = similar(P)

    init_gs_deep_quench!(dim, Q, P, ρ)

    H_hop = hopping_Hamiltonian_1D_dense(dim)

    subfolder = "$(save_base)run_$(seed)/"
    mkpath(subfolder)

    for step in 1:nS
        time_ms = @elapsed begin
            evolve!(
                dim,
                Q,
                P,
                ρ,
                dt,
                r,
                H_hop,
                λ,
                ρ2,
                Kρ_sum,
                Q2,
                KQ_sum,
                P2,
                KP_sum
            )
        end * 1000

        println("Step $step time (ms): ", round(time_ms, digits=2))

        if step % 100 == 0
            println("== Step $step Report ==")
            println("  Elapsed time per step: $(round(time_ms, digits=2)) ms")
            println("  GPU memory usage:")
            CUDA.memory_status()
        end

        if step % save_interval == 0
            save_vector_line(subfolder * "Q.txt", Array(Q))
            save_vector_line(subfolder * "P.txt", Array(P))
            save_vector_line(subfolder * "n.txt", real.(Array(diag(ρ))))
        end
    end
end


# ============================================================
# Run experiment from SLURM job id
# ============================================================

function run_experiments()
    slurm_job_id = get(ENV, "SLURM_JOB_ID", "default")

    if slurm_job_id == "default"
        seed = 1
    else
        seed = parse(Int, slurm_job_id)
    end

    save_base = "/standard/mott-physics/Jason/transformer/Holstein/force_prediction/data/r0.300_lambda0.600_dim1024_dt0.01_nS1200000/"

    @time main(seed, save_base)
end


run_experiments()
