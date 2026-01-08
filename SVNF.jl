################################################################################
# Allan Eduardo Feitosa
# 27 Jun 2025
# SVNF Algorithm

function
    SVNF!( 
        bitvector::Vector{Bool},
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},
        Lf::Vector{Float64},
        Nc::Vector{Vector{Int}},
        Nv::Vector{Vector{Int}},
        phi::Union{Vector{Float64},Nothing},
        msum_factor::Union{Float64,Nothing},
        msum2::Bool,
        num_reps::Int,
        newLr::Matrix{Float64},
        Residues::Matrix{Float64},
        rbp_not_converged::Bool,
        twoLS::Int,
        N::Int      
    )

    count = 0
    count_zeros = 0

    @inbounds while count < num_reps

        # display("count = $count")
    
        for vn in 0:N-1

            # count += 1

            vj = rem(vn + twoLS,N) + 1   # jump to non_punctured nodes

            # display("vj = $vj")

            # 1) Find largest residue  and coordenates
            Nvj = Nv[vj]
            cimax, vjmax = findmaxedge_SVNF(Residues,vj,Nvj,Nc)
            if cimax ≠ 0
                count_zeros = 0

                # 2) update check to node message Lr[cimax,vjmax]
                Ncimax = Nc[cimax]
                limax = LinearIndices(Lr)[cimax,vjmax]
                if msum2
                    Lr[limax] = calc_Lr_no_opt(Ncimax,cimax,vjmax,Lq)
                else
                    Lr[limax] = newLr[limax]
                end
                
                count += 1
                
                # 3) set maximum residue to zero
                Residues[limax] = 0.0

                # 4) Calculate Ld of vjmax and bitvector[vjmax]
                Nvjmax = Nv[vjmax]
                Ld = calc_Ld(vjmax,Nvjmax,Lf,Lr)
                bitvector[vjmax] = signbit(Ld)

                for ci in Nvjmax
                    if ci ≠ cimax
                        # 6) update Nv messages Lq[ci,vjmax]
                        li = LinearIndices(Lq)[ci,vjmax]
                        Lq[li] = tanhLq(Ld,Lr[li],msum_factor)
                        # 7) calculate Residues
                        Nci = Nc[ci]    
                        for vj in Nci
                            if vj ≠ vjmax
                                li = LinearIndices(Lr)[ci,vj]
                                newlr = calc_Lr(Nci,ci,vj,Lq,msum_factor)
                                newLr[li] = newlr
                                Residues[li] = abs(newlr - Lr[li])
                            end
                        end
                    end
                end
            else
                count_zeros += 1        
            end
            if count > num_reps
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