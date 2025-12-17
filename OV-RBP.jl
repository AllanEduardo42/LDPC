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
        num_reps::Int,
        Residues::Vector{Float64},
        rbp_not_converged::Bool,
        msum_factor::Union{Float64,Nothing},
        newLv::Vector{Float64},
        Lv::Vector{Float64},
        C::Vector{Bool},
        upc::Vector{Int},
        max_upc::Int,
        size_C::Int
    )   

    @fastmath @inbounds for e in 1:num_reps

        # display("e = $e")

        # Dynamic selection

        vjmax = 0
        maxresidue = 0.0
        
        if sum(C) > 0

            for vj in eachindex(upc)
                if upc[vj] == max_upc
                    if C[vj]
                        residue = Residues[vj] 
                        if residue > maxresidue
                            maxresidue = residue
                            vjmax = vj
                        end
                    end
                end
            end

            if vjmax == 0 # C i
                max_upc = 0
                for p in upc
                    if p > max_upc
                        max_upc = p
                    end
                end
                if max_upc == 0
                    max_upc = 1
                end    
                for vj in eachindex(Nv)
                    if C[vj]
                        residue = Residues[vj]
                        if residue > maxresidue
                            maxresidue = residue
                            vjmax = vj 
                        end
                    end
                end  
            end
        else
            for vj in eachindex(Nv)
                residue = Residues[vj]
                if residue > maxresidue
                    maxresidue = residue
                    vjmax = vj 
                end
            end
        end

        if vjmax == 0
            # # display("bp")
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
            newlvOsc = oldlv + newLv[vjmax]
            # if sign(oldlv)*sign(newlvOsc) < 0 
            #     C[vjmax] = true
            # end
            Lv[vjmax] = newlvOsc
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
                    end
                    upc[vj] = 0
                    li = LinearIndices(Lr)[ci,vj]
                    newLr[li] = calc_Lr(Nci,ci,vj,Lq,msum_factor)
                    oldlv = Lv[vj]
                    newlv = calc_Ld(vj,Nv[vj],Lf,newLr)
                    Residues[vj] = abs(newlv - oldlv)        
                    newLv[vj] = newlv
                    if sign(oldlv)*sign(newlv) < 0
                        C[vj] = true
                        size_C += 1
                        for cii in Nv[vj]
                            if _calc_syndrome(bitvector,Nc[cii])
                                upc[vj] += 1
                            end
                        end
                    end
                end
            end
        end
    end

    return rbp_not_converged, max_upc, size_C

end










        





