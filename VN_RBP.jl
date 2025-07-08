################################################################################
# Allan Eduardo Feitosa
# 31 Mar 2025
# RBP Sum-Product Algorithm with residual decaying factor

include("./RBP functions/findmaxnode.jl")
include("./RBP functions/calc_residue.jl")

function
    VN_RBP!(
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
        Factors::Vector{Float64},
        alpha::Vector{Float64},
        rbp_not_converged::Bool,
        raw::Bool
    )

    @fastmath @inbounds if raw
        for e in 1:num_reps

            # display("e = $e")

            vjmax = findmaxnode(alpha)
            if vjmax == 0
                rbp_not_converged = false
                break # i.e., BP has converged
            end
            Factors[vjmax] *= decayfactor
            alpha[vjmax] = 0.0

            Nvjmax = Nv[vjmax]
            for ci in Nvjmax
                li = LinearIndices(Lr)[ci,vjmax]
                Lr[li] = newLr[li]
            end 

            for ci in Nvjmax
                # 5) update Nv messages Lq[ci,vnmax]
                Lq[ci,vjmax] = calc_Lq(Nvjmax,ci,vjmax,Lr,Lf)
                # 6) calculate alpha
                Nci = Nc[ci]
                for vj in Nci
                    if vj ≠ vjmax
                        newlr = calc_Lr(Nci,ci,vj,Lq)
                        li = LinearIndices(Lr)[ci,vj]
                        newLr[li] = newlr
                        alpha[vj] = calc_residue_VN_NW_raw(newlr,Lr[li],Factors[vj])
                    end
                end
            end
        end

        # 7) update bitvector
        for vj in eachindex(Nv)
            Ld = calc_Ld(vj,Nv[vj],Lf,Lr)
            bitvector[vj] = signbit(Ld)
        end
    
    else

        for e in 1:num_reps

            # display("e = $e")

            vjmax = findmaxnode(alpha)
            if vjmax == 0
                rbp_not_converged = false
                break # i.e., BP has converged
            end
            Factors[vjmax] *= decayfactor
            alpha[vjmax] = 0.0

            Nvjmax = Nv[vjmax]
            for ci in Nvjmax
                li = LinearIndices(Lr)[ci,vjmax]
                RBP_update_Lr!(li,Lr,newLr,ci,vjmax,Nc[ci],Lq,signs,phi)
            end   
            
            Ld = calc_Ld(vjmax,Nvjmax,Lf,Lr)
            bitvector[vjmax] = signbit(Ld)

            for ci in Nvjmax
                # update Nv messages Lq[ci,vjmax]
                li = LinearIndices(Lq)[ci,vjmax]
                Lq[li] = tanhLq(Ld - Lr[li],signs)
                # calculate the new check to node messages
                Nci = Nc[ci]
                A, B, C, D = calc_ABCD!(Lq,ci,Nci,signs,phi)
                for vj in Nci
                    if vj ≠ vjmax
                        li = LinearIndices(newLr)[ci,vj]
                        newlr = calc_Lr(A,B,C,D,vj,Lq[li],signs,phi)
                        newLr[li] = newlr
                        residue = abs(newlr - Lr[li])*Factors[vj]
                        if residue > alpha[vj]
                            alpha[vj] = residue
                        end
                    end
                end
            end
        end
    end

    return rbp_not_converged

end
