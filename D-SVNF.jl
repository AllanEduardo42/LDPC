################################################################################
# Allan Eduardo Feitosa
# 6 Jan 2026
# D-SVNF Algorithm

function
    D_SVNF!( 
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
        N::Int,        
        u::Vector{Bool}
    )

    count = 0
    vp = 1
    u[1] = true

    @fastmath @inbounds while count < num_reps

        # display("e = $count")

        count += 1

        Nvp = Nv[vp]
        # display("vp = $vp")
        # display("Nvp = $Nvp")
        maxresidue = 0.0
        cimax = 0
        vjmax = 0
        no_unsat_check = true
        for ci in Nvp
            Nci = Nc[ci] 
            if _calc_syndrome(bitvector,Nci)
                # display("s[$ci] = true")  
                if no_unsat_check
                    no_unsat_check = false
                end   
                for vj in Nci
                    if vj ≠ vp
                        # display("   vj = $vj")
                        li = LinearIndices(Lr)[ci,vj]
                        oldlr = Lr[li]
                        newlr = calc_Lr(Nci,ci,vj,Lq,msum_factor)
                        newLr[li] = newlr
                        residue = abs(newlr - oldlr)
                        if residue ≥ maxresidue
                            maxresidue = residue
                            cimax = ci
                            vjmax = vj
                        end
                        # display("   residue = $residue")
                    end
                end
            end
        end            
        if no_unsat_check
            # display("no unsat check")
            for ci in Nvp
                Nci = Nc[ci] 
                for vj in Nci
                    if vj ≠ vp
                        # display("   vj = $vj")
                        li = LinearIndices(Lr)[ci,vj]
                        oldlr = Lr[li]
                        newlr = calc_Lr(Nci,ci,vj,Lq,msum_factor)
                        newLr[li] = newlr
                        residue = abs(newlr - oldlr)
                        if residue ≥ maxresidue
                            maxresidue = residue
                            cimax = ci
                            vjmax = vj
                        end
                        # display("   residue = $residue")
                    end
                end
            end 
        end

        # display("maxresidue = $maxresidue")
        # display("cimax = $cimax")
        # display("vjmax = $vjmax")


        Ncimax = Nc[cimax]
        limax = LinearIndices(Lr)[cimax,vjmax]
        if msum2
            Lr[limax] = calc_Lr_no_opt(Ncimax,cimax,vjmax,Lq)
        else
            Lr[limax] = newLr[limax]
        end
        Residues[limax] = 0.0

        Nvjmax = Nv[vjmax]
        Ld = calc_Ld(vjmax,Nvjmax,Lf,Lr)
        bitvector[vjmax] = signbit(Ld)

        for ci in Nvjmax
            if ci ≠ cimax
                # 6) update V2C messages Lq[ci,vjmax]
                li = LinearIndices(Lq)[ci,vjmax]
                Lq[li] = tanhLq(Ld,Lr[li],msum_factor)
            end
        end

        # display("u[vjmax] = $(u[vjmax])")
        if u[vjmax] 
            reset_u = true
            for outer vp in eachindex(u)
                if !u[vp]
                    reset_u = false
                    break
                end
            end

            # display("   vp = $vp")
            # display("   reset = $reset_u")

            if reset_u
                u .= false
                vp = 1
            end
        else
            vp = vjmax
        end

        u[vp] = true

        # println()
        # println()
    end

    return rbp_not_converged

end