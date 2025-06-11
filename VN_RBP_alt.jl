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
        num_edges::Int,
        newLr::Matrix{Float64},
        Factors::Matrix{Float64},
        coords::Matrix{Int},
        indices::Matrix{Int},
        residues::Vector{Float64},
        relative::Bool,
        rbp_not_converged::Bool,
        alpha::Vector{Float64},
        ci_alpha::Vector{Int}
    )

    @fastmath @inbounds for e in 1:num_edges

        # display("### e = $e")
        vjmax = findmaxnode(alpha)
        if vjmax == 0
            rbp_not_converged = false
            break # i.e., BP has converged
        end

        # 2) Decay the RBP factor corresponding to the maximum residue
        # Factors[limax] *= decayfactor

        # 3) update check to node message Lr[cimax,vjmax]
        Nvjmax = Nv[vjmax]
        for ci in Nvjmax
            li = LinearIndices(Lr)[ci,vjmax]
            RBP_update_Lr!(li,Lr,newLr,ci,vjmax,Nc[ci],Lq,aux,signs,phi)
            # 4) set maximum residue to zero
            residues[indices[li]] = 0.0
        end
        alpha[vjmax] = 0.0

        # 5) Calculate Ld of vjmax and bitvector[vjmax]
        Ld = calc_Ld(vjmax,Nvjmax,Lf,Lr)
        bitvector[vjmax] = signbit(Ld)

        for ci in Nvjmax
            # if ci ≠ cimax
                # 6) update Nv messages Lq[ci,vjmax]
                li = LinearIndices(Lq)[ci,vjmax]
                Lq[li] = Ld - Lr[li]
                # 7) calculate residues
                Nci = Nc[ci]
                # display("Nci = $Nci")    
                A, B, C, D = calc_ABCD!(aux,signs,phi,Lq,ci,Nci)
                for vj in Nci
                    # display("vj = $vj")
                    if vj ≠ vjmax
                        newlr = calc_Lr(A,B,C,D,vj,aux,signs,phi)
                        li = LinearIndices(Lr)[ci,vj]
                        newLr[li] = newlr                                              
                        residue = calc_residue(newlr,Lr[li],Factors[li],
                                                        relative,Lq[li])
                        residues[indices[li]] = residue
                        if ci == ci_alpha[vj]
                            if residue > alpha[vj]
                                alpha[vj] = residue
                            else
                                alp = residue
                                ci_alp = ci
                                for ca in Nv[vj]
                                    if ca ≠ ci
                                        residue = residues[indices[ca,vj]]
                                        if residue > alp
                                            alp = residue
                                            ci_alp = ca
                                        end
                                    end
                                end
                                alpha[vj]= alp
                                ci_alpha[vj] = ci_alp
                            end
                        elseif residue > alpha[vj]
                            alpha[vj]= residue
                            ci_alpha[vj] = ci
                        end
                    end
                end
            # end
        end
    end

    return rbp_not_converged
end

