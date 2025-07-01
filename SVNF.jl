################################################################################
# Allan Eduardo Feitosa
# 27 Jun 2025
# SVNF Algorithm with residual decaying factor

function
    SVNF!( 
        bitvector::Vector{Bool},
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},
        Lf::Vector{Float64},
        Nc::Vector{Vector{Int}},
        Nv::Vector{Vector{Int}},
        aux::Vector{Float64},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{Float64},Nothing},
        N::Int,
        num_edges::Int,
        newLr::Matrix{Float64},
        residues::Matrix{Float64},
        rbp_not_converged::Bool,
        twoLS::Int
        )

    count = 0
    count_zeros = 0
    @fastmath @inbounds while count < num_edges
    
        for n in 0:N-1

            vj = rem(n + twoLS,N) + 1

            # display("vj = $vj")

            # 1) Find largest residue  and coordenates
            Nvj = Nv[vj]
            cimax, vjmax = findmaxedge_SVNF(residues,vj,Nvj,Nc)
            if cimax ≠ 0
                count_zeros = 0
                rbp_not_converged = true
                #     rbp_not_converged = false
                #     break # i.e., BP has converged
                # end

                # 2) Decay the RBP factor corresponding to the maximum residue
                # Factors[limax] *= decayfactor

                # 3) update check to node message Lr[cimax,vjmax]
                Ncimax = Nc[cimax]
                limax = LinearIndices(Lr)[cimax,vjmax]
                RBP_update_Lr!(limax,Lr,newLr,cimax,vjmax,Ncimax,Lq,aux,signs,phi)
                count += 1

                # 4) set maximum residue to zero
                residues[limax] = 0.0

                # 5) Calculate Ld of vjmax and bitvector[vjmax]
                Nvjmax = Nv[vjmax]
                Ld = calc_Ld(vjmax,Nvjmax,Lf,Lr)
                bitvector[vjmax] = signbit(Ld)

                for ci in Nvjmax
                    if ci ≠ cimax
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
                                residues[li] = abs(newlr - Lr[li])
                            end
                        end
                    end
                end
            else
                count_zeros += 1                
            end
            if count > num_edges
                break
            end
        end
        if count_zeros ≥ N
            rbp_not_converged = false
            break
        end
    end

    return rbp_not_converged

end