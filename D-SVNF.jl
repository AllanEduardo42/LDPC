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
        N::Int,        
        u::Vector{Bool}
    )

    count = 1
    count_zeros = 0
    vp = 1

    @fastmath @inbounds while count <= num_reps

        u[vp] = true

        # display("e = $count")

        Nvp = Nv[vp]
        # display("   vp = $vp")
        # display("   Nvp = $Nvp")
        maxresidue = 0.0
        cimax = 0
        vjmax = 0
        no_unsat_check = true
        for ci in Nvp
            Nci = Nc[ci] 
            if _calc_syndrome(bitvector,Nci)
                # display("       s[$ci] = true")  
                if no_unsat_check
                    no_unsat_check = false
                end   
                maxresidue, cimax, vjmax = calc_newLr!(Lq,Lr,ci,Nci,msum_factor,newLr,
                                                maxresidue,cimax,vjmax,vp)
            end
        end            
        if no_unsat_check
            # display("       no unsatisfied check") 
            for ci in Nvp
                # display("       ci = $ci")
                Nci = Nc[ci] 
                maxresidue, cimax, vjmax = calc_newLr!(Lq,Lr,ci,Nci,msum_factor,newLr,
                                                maxresidue,cimax,vjmax,vp)
            end 
        end

        # display("   maxresidue = $maxresidue")
        # display("   cimax = $cimax")
        # display("   vjmax = $vjmax")

        if vjmax == 0
            vp = find_vp!(u)  
            count_zeros += 1      
        else
            count_zeros = 0
            count += 1
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

            # display("   u[vjmax] = $(u[vjmax])")
            if u[vjmax] 
                vp = find_vp!(u)
            else
                vp = vjmax
            end
        end
        if count_zeros ≥ N
            rbp_not_converged = false
            break
        end    
    end

    return rbp_not_converged

end

function 
    find_vp!(
        u::Vector{Bool}
    )

    # print_test("u",u)

    reset_u = true
    vp = 1
    for outer vp in eachindex(u)
        if !u[vp]
            reset_u = false
            break
        end
    end

    # display("       vp = $vp")
    # display("       reset = $reset_u")

    if reset_u
        u .= false
        vp = 1
    end

    return vp

end

function 
    calc_newLr!(
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},
        ci::Int,
        Nci::Vector{Int},
        msum_factor::Union{Float64,Nothing},
        newLr::Matrix{Float64},
        maxresidue::Float64,
        cimax::Int,
        vjmax::Int,
        vp::Int
    )
    
    for vj in Nci
        if vj ≠ vp
            # display("           vj = $vj")
            li = LinearIndices(Lr)[ci,vj]
            newlr = calc_Lr(Nci,ci,vj,Lq,msum_factor)
            newLr[li] = newlr
            residue = abs(newlr - Lr[li])
            if residue > maxresidue
                maxresidue = residue
                cimax = ci
                vjmax = vj
            end
            # display("           residue = $residue")
        end
    end

    return maxresidue, cimax, vjmax
end