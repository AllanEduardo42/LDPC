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

        # println()
        # display("### e = $e")

        # 1) Find largest residue  and coordenates
        # max_edge = findmaxedge(residues)
        # if max_edge == 0.0
        #     rbp_not_converged = false
        #     break # i.e., BP has converged
        # else
        #     cimax = coords[1,max_edge]
        #     vjmax = coords[2,max_edge]
        #     limax = coords[3,max_edge]
        # end
        vjmax = findmaxnode(alpha)
        if vjmax == 0
            rbp_not_converged = false
            break # i.e., BP has converged
        end

        # display("(cimax,vjmax) = ($cimax,$vjmax)")

        # 2) Decay the RBP factor corresponding to the maximum residue
        # Factors[limax] *= decayfactor

        # 3) update check to node message Lr[cimax,vjmax]
        Nvjmax = Nv[vjmax]
        # display("Nvjmax = $Nvjmax")
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
            # println()
            # display("ci = $ci")
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
                        # display("ci_alpha[$vj] = $(ci_alpha[vj])")
                        if ci == ci_alpha[vj]
                            # display("residue = $residue")
                            # display("alpha[$vj] = $(alpha[vj])")
                            if residue > alpha[vj]
                                alpha[vj] = residue
                            else
                                alp = residue
                                ci_alp = ci
                                for ca in Nv[vj]
                                    if ca ≠ ci
                                        residue = residues[indices[ca,vj]]
                                        # display("           new residue = $residue")
                                        if residue > alp
                                            alp = residue
                                            ci_alp = ca
                                        end
                                    end
                                end
                                alpha[vj]= alp
                                ci_alpha[vj] = ci_alp
                            end
                            # display("alpha[$vj] = $(alpha[vj])")
                            # display("ci_alpha[$vj] = $(ci_alpha[vj])")
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

# RAW
function
    VN_RBP_RBP!(
        bitvector::Vector{Bool},
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},
        Lf::Vector{Float64},
        Nc::Vector{Vector{Int}},
        Nv::Vector{Vector{Int}},
        ::Nothing,
        ::Nothing,
        ::Nothing,
        decayfactor::Float64,
        num_edges::Int,
        newLr::Matrix{Float64},
        Factors::Matrix{Float64},
        coords::Matrix{Int},
        indices::Matrix{<:Integer},
        residues::Vector{Float64},
        relative::Bool,
        rbp_not_converged::Bool
    )

    @fastmath @inbounds for e in 1:num_edges

        # # display("e = $e")

        # 1) Find largest residue  and coordenates
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
        # 3) update check to node message Lr[cnmax,vnmax]
        Lr[limax] = newLr[limax]
        # 4) set maximum residue to zero
        residues[max_edge] = 0.0

        Nvjmax = Nv[vjmax]
        for ci in Nvjmax
            if ci ≠ cimax
                # 5) update Nv messages Lq[ci,vnmax]
                Lq[ci,vjmax] = calc_Lq(Nvjmax,ci,vjmax,Lr,Lf)
                # 6) calculate residues
                Nci = Nc[ci]
                for vj in Nci
                    if vj ≠ vjmax
                        newlr = calc_Lr(Nci,ci,vj,Lq)
                        li = LinearIndices(Lr)[ci,vj]
                        newLr[li] = newlr
                        residue = calc_residue_raw(newlr,Lr[li],Factors[li],
                                                   relative,Lq[li])
                        residues[indices[li]] = residue
                    end
                end
            end
        end
    end

    # 7) update bitvector
    for vj in eachindex(Nv)
        ci = Nv[vj][1]
        li = LinearIndices(Lr)[ci,vj]
        bitvector[vj] = signbit(Lr[li] + Lq[li])
    end

    return rbp_not_converged
end

