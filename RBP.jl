################################################################################
# Allan Eduardo Feitosa
# 07 Apr 2025
# RBP and List-RBP Algorithm with residual decaying factor

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/findmaxedge.jl")
include("./RBP functions/calc_residue.jl")

function
    RBP!(
        bitvector::Vector{Bool},
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},
        Lf::Vector{Float64},
        Nc::Vector{Vector{Int}},
        Nv::Vector{Vector{Int}},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{Float64},Nothing},
        decayfactor::Float64,
        num_reps::Int,
        newLr::Matrix{Float64},
        alpha::Vector{Float64},
        Residues::Matrix{Float64},
        Factors::Matrix{Float64},
        rbp_not_converged::Bool,
        two_ways::Bool,
        correction::Bool,
        switch_C_VN::Bool
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
        if correction
            for ci in Nvjmax
                li = LinearIndices(Lr)[ci,vjmax]
                # 2) Decay the RBP factor corresponding to the maximum residue
                Factors[li] *= decayfactor
                # 3) update check to node message Lr[cnmax,vnmax]
                Lr[li] = newLr[li]
                # 4) set maximum residue to zero
                Residues[li] = 0.0
            end
        else
            limax = LinearIndices(Lr)[cimax,vjmax]
            # 2) Decay the RBP factor corresponding to the maximum residue
            Factors[limax] *= decayfactor
            # 3) update check to node message Lr[cnmax,vnmax]
            Lr[limax] = newLr[limax]
            # 4) set maximum residue to zero
            Residues[limax] = 0.0
        end       

        Ld = calc_Ld(vjmax,Nvjmax,Lf,Lr)
        bitvector[vjmax] = signbit(Ld)

        for ci in Nvjmax
            if switch_C_VN || ci ≠ cimax
                # 5) update Nv messages Lq[ci,vnmax]
                li = LinearIndices(Lq)[ci,vjmax]
                alp = Residues[li]
                Lq[li] = tanh(0.5*(Ld - Lr[li]))
                # 6) calculate Residues
                Nci = Nc[ci]
                for vj in Nci
                    if vj ≠ vjmax
                        li = LinearIndices(Lr)[ci,vj]
                        alp, residue = calc_residue!(Lq,Lr,newLr,li,ci,vj,Nci,Factors[li],alp)
                        Residues[li] = residue
                    end
                end
                alpha[ci] = alp
            end
        end

        if two_ways
            if !correction
                for ci in Nvjmax
                    if ci ≠ cimax
                        li = LinearIndices(Lr)[ci,vjmax]
                        Lr[li] = newLr[li]
                    end
                end
            end
            alp = 0.0
            # 5) update Nv messages Lq[cimax,vnmax]
            li = LinearIndices(Lq)[cimax,vjmax]
            Lq[li] = tanh(0.5*(Ld - Lr[li]))
            # 6) calculate Residues
            Ncimax = Nc[cimax]
            for vj in Ncimax
                if vj ≠ vjmax
                    li = LinearIndices(Lr)[cimax,vj]
                    alp, residue = calc_residue!(Lq,Lr,newLr,li,cimax,vj,Ncimax,Factors[li],alp)
                    Residues[li] = residue
                end
            end
            alpha[cimax] = alp
        end
    end

    return rbp_not_converged
end

