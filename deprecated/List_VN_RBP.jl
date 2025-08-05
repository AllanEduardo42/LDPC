################################################################################
# Allan Eduardo Feitosa
# 31 Mar 2025
# RBP Sum-Product Algorithm with residual decaying factor

include("./RBP functions/findmaxnode.jl")
include("./RBP functions/calc_residue.jl")

function
    List_VN_RBP!(
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
        coords::Vector{Int},
        inlist::Vector{Bool},
        listsize::Int,
        rbp_not_converged::Bool
    )

    @fastmath @inbounds for e in 1:num_reps

        if alpha[1] == 0.0
            init_list_VN_RBP!(Lq,Nc,Nv,newLr,Factors,inlist,
                                                        alpha,coords,listsize)
            if alpha[1] == 0.0
                rbp_not_converged = false
                break
            end
        end
        vjmax = coords[1]
        inlist[vjmax] = false

        # update the list
        for i in 1:listsize
            alpha[i] = alpha[i+1]
            coords[i] = coords[i+1]
        end 

        Factors[vjmax] *= decayfactor        

        Nvjmax = Nv[vjmax]
        for ci in Nvjmax
            li = LinearIndices(Lr)[ci,vjmax]
            Lr[li] = newLr[li]
        end   
        
        Ld = calc_Ld(vjmax,Nvjmax,Lf,Lr)
        bitvector[vjmax] = signbit(Ld)

        for ci in Nvjmax
            # update Nv messages Lq[ci,vjmax]
            li = LinearIndices(Lq)[ci,vjmax]
            Lq[li] = tanh(0.5*(Ld - Lr[li]))
            # calculate the new check to node messages
            Nci = Nc[ci]
            for vj in Nci
                if vj ≠ vjmax
                    li = LinearIndices(Lr)[ci,vj]
                    newlr = calc_Lr(Nci,ci,vj,Lq)
                    newLr[li] = newlr
                    residue = abs(newlr - Lr[li])*Factors[vj]
                    if inlist[vj]
                        pos = 0
                        for i = 1:listsize
                            if coords[i] == vj
                                pos = i
                                break
                            end
                        end
                        if pos == 0
                            throw(error("($(coords[1,i]),$(coords[1,i])) is registered as being on the list, but it's not."))
                        end
                        if residue > alpha[pos]
                            # remove from list
                            inlist[vj] = false
                            # update the list
                            for i in pos:listsize
                                alpha[i] = alpha[i+1]
                                coords[i] = coords[i+1]
                            end
                            add_list_VN!(residue,alpha,listsize,inlist,coords,vj)
                        end
                    else
                        add_list_VN!(residue,alpha,listsize,inlist,coords,vj)
                    end                     
                end
            end
        end
    end

    return rbp_not_converged

end

function 
    add_list_VN!(
        residue::Float64,
        alpha::Vector{Float64},
        listsize::Int,
        inlist::Vector{Bool},
        coords::Vector{Int},
        vj::Int
    )

    if residue > alpha[listsize]
        if residue ≥ alpha[1]
            i = 1
        else
            d = listsize >> 1
            i = d
            while d > 1
                d >>= 1
                if residue ≥ alpha[i]
                    i -= d
                else
                    i += d
                end
            end
            if residue < alpha[i]
                i += 1
            end
        end

        last = coords[end-1]
        if last ≠ 0
            inlist[last] = false
        end
        inlist[vj] = true

        for j=listsize:-1:i+1
            alpha[j] = alpha[j-1]
            coords[j] = coords[j-1]
        end
        coords[i] = vj
        alpha[i] = residue        
    end
end
