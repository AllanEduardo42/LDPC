function
    OV_RBP!(
        bitvector::Vector{Bool},
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},
        newLr::Matrix{Float64},
        Lf::Vector{Float64},
        Nc::Vector{Vector{Int}},
        Nv::Vector{Vector{Int}},
        phi::Union{Vector{Float64},Nothing},
        msum_factor::Union{Float64,Nothing},
        msum2::Bool,
        num_reps::Int,
        newLv::Vector{Float64},
        Residues::Vector{Float64},
        rbp_not_converged::Bool,
        Lv::Vector{Float64},
        C::Vector{Bool},
        upc::Vector{Int},
        max_upc::Int,
        size_C::Int
    )   

    @fastmath @inbounds for e in 1:num_reps

        # println()
        # display("e = $e")
        # display("P = $max_upc")

        # Dynamic selection

        vjmax = 0
        maxresidue = 0.0
        
        if size_C > 0 # if C is not empty

            for vj in eachindex(upc)    # find vjmax in C1, if C1 is not empty
                if upc[vj] == max_upc
                    if C[vj]
                        # print("$vj ")
                        residue = Residues[vj] 
                        # println("residue = $residue")
                        if residue > maxresidue
                            maxresidue = residue
                            vjmax = vj
                        end
                    end
                end
            end
            # println()

            if vjmax == 0 # if C1 is empty, find vjmax in C2   
                # display("C2")
                for vj in eachindex(Nv)
                    if C[vj]
                        residue = Residues[vj]
                        if residue > maxresidue
                            maxresidue = residue
                            vjmax = vj 
                        end
                    end
                end 
                # Find new maximum number of unsatisfied paity check equations
                max_upc = 1
                for p in upc
                    if p > max_upc
                        max_upc = p
                    end
                end
            end
        else # if C is empty
            for vj in eachindex(Nv)
                residue = Residues[vj]
                if residue > maxresidue
                    maxresidue = residue
                    vjmax = vj 
                end
            end
        end

        # display("vjmax = $vjmax")
        # display("maxresidue = $maxresidue")

        if vjmax == 0
            rbp_not_converged = false
            break # i.e., BP has converged
        end

        Nvjmax = Nv[vjmax]

        # C2V

        for ci in Nvjmax
            li = LinearIndices(Lr)[ci,vjmax]
            Lr[li] = newLr[li]
        end

        # Oscillation decision
        if C[vjmax]
            C[vjmax] = false
            size_C -= 1
            upc[vjmax] = 0
            oldlv = Lv[vjmax]
            Lv[vjmax] = oldlv + newLv[vjmax]
            # Lv[vjmax] = newLv[vjmax]
        else
            Lv[vjmax] = newLv[vjmax]
        end

        bitvector[vjmax] = signbit(Lv[vjmax])

        Residues[vjmax] = 0.0

        # V2C

        for ci in Nvjmax
            li = LinearIndices(Lr)[ci,vjmax]
            Lq[li] = tanhLq(Lv[vjmax] - Lr[li],msum_factor)
        end

        # Lv

        for ci in Nvjmax
            Nci = Nc[ci]
            for vj in Nci
                if vj != vjmax                    
                    if C[vj]
                        C[vj] = false
                        size_C -= 1
                        upc[vj] = 0
                    end
                    li = LinearIndices(Lr)[ci,vj]
                    newLr[li] = calc_Lr(Nci,ci,vj,Lq,msum_factor)
                    oldlv = Lv[vj]
                    newlv = calc_Ld(vj,Nv[vj],Lf,newLr)
                    Residues[vj] = abs(newlv - oldlv)        
                    newLv[vj] = newlv
                    if sign(oldlv)*sign(newlv) < 0
                        C[vj] = true
                        size_C += 1
                        count = 0
                        for cii in Nv[vj]
                            if _calc_syndrome(bitvector,Nc[cii])
                                count += 1
                            end
                        end
                        upc[vj] = count
                    end
                end
            end
        end
    end

    return rbp_not_converged, max_upc, size_C

end










        





