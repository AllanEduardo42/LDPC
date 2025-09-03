################################################################################
# Allan Eduardo Feitosa
# 07 Apr 2025
# RBP and List-RBP Algorithm with residual decaying factor

include("./RBP functions/findmaxedge.jl")
include("./RBP functions/calc_residue.jl")

function
    C_R_RBP!(
        bitvector::Vector{Bool},
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},
        Lf::Vector{Float64},
        Nc::Vector{Vector{Int}},
        Nv::Vector{Vector{Int}},
        phi::Union{Vector{Float64},Nothing},
        decayfactor::Float64,
        num_reps::Int,
        newLr::Matrix{Float64},
        alpha::Vector{Float64},
        Residues::Matrix{Float64},
        Factors::Matrix{Float64},
        rbp_not_converged::Bool,
        msum_factor::Union{Float64,Nothing},
        msum2::Bool
    )
    
    @fastmath @inbounds for e in 1:num_reps

        # display("e = $e")

        # 1) Find largest residue  and coordenates
        cimax, vjmax = findmaxedge(Residues,alpha,Nc)
        if cimax == 0.0
            rbp_not_converged = false
            break # i.e., BP has converged
        end

        Nvjmax = Nv[vjmax]
        for ci in Nvjmax
            li = LinearIndices(Lr)[ci,vjmax]
            # 2) Decay the RBP factor corresponding to the maximum residue
            Factors[li] *= decayfactor
            # 3) update check to node message Lr[cnmax,vnmax]
            if msum2
                Lr[li] = calc_Lr_no_opt(Nc[ci],ci,vjmax,Lq)
            else
                Lr[li] = newLr[li]
            end
            # 4) set maximum residue to zero
            Residues[li] = 0.0
        end

        Ld = calc_Ld(vjmax,Nvjmax,Lf,Lr)
        bitvector[vjmax] = signbit(Ld)

        for ci in Nvjmax
            # 5) update Nv messages Lq[ci,vnmax]
            li = LinearIndices(Lq)[ci,vjmax]
            alp = Residues[li]
            Lq[li] = tanhLq(Ld - Lr[li],msum_factor)
            # 6) calculate Residues
            Nci = Nc[ci]
            for vj in Nci
                if vj ≠ vjmax
                    li = LinearIndices(Lr)[ci,vj]
                    alp, residue = calc_residue!(Lq,Lr,newLr,li,ci,vj,Nci,Factors[li],alp,msum_factor)
                    Residues[li] = residue
                end
            end
            alpha[ci] = alp
        end
    end

    return rbp_not_converged
end