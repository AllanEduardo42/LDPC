################################################################################
# Allan Eduardo Feitosa
# 07 Apr 2025
# RBP and List-RBP Algorithm with residual decaying factor

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/findmaxedge.jl")
include("./RBP functions/calc_residue.jl")

function
    VN_RBP_ALT!(
        bitvector::Vector{Bool},
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},
        Lf::Vector{Float64},
        Nc::Vector{Vector{Int}},
        Nv::Vector{Vector{Int}},
        aux::Vector{Float64},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{Float64},Nothing},
        decayfactor::Float64,
        num_reps::Int,
        newLr::Matrix{Float64},
        Factors::Matrix{Float64},
        coords::Matrix{Int},
        indices::Matrix{Int},
        residues::Vector{Float64},
        rbp_not_converged::Bool
    )

    @fastmath @inbounds for e in 1:num_reps

        # display("### e = $e")
        max_edge = findmaxedge(residues)
        if max_edge == 0.0
            rbp_not_converged = false
            break # i.e., BP has converged
        else
            cimax = coords[1,max_edge]
            vjmax = coords[2,max_edge]
            limax = coords[3,max_edge]
        end

        # 2) Decay the RBP factor corresponding to the maximum residue
        Factors[limax] *= decayfactor

        # 3) update check to node message Lr[cimax,vjmax]
        Nvjmax = Nv[vjmax]
        for ci in Nvjmax
            li = LinearIndices(Lr)[ci,vjmax]
            RBP_update_Lr!(li,Lr,newLr,ci,vjmax,Nc[ci],Lq,aux,signs,phi)
            # 4) set maximum residue to zero
            residues[indices[li]] = 0.0
        end

        # 5) Calculate Ld of vjmax and bitvector[vjmax]
        Ld = calc_Ld(vjmax,Nvjmax,Lf,Lr)
        bitvector[vjmax] = signbit(Ld)

        for ci in Nvjmax
            # 6) update Nv messages Lq[ci,vjmax]
            li = LinearIndices(Lq)[ci,vjmax]
            Lq[li] = Ld - Lr[li]
            # 7) calculate residues
            Nci = Nc[ci] 
            A, B, C, D = calc_ABCD!(aux,signs,phi,Lq,ci,Nci)
            for vj in Nci
                if vj ≠ vjmax
                    newlr = calc_Lr(A,B,C,D,vj,aux,signs,phi)
                    li = LinearIndices(Lr)[ci,vj]
                    newLr[li] = newlr                                              
                    residue = calc_residue(newlr,Lr[li],Factors[li],
                                                    false,Lq[li])
                    residues[indices[li]] = residue*Factors[li]
                end
            end
        end
    end

    return rbp_not_converged
end

