################################################################################
# Allan Eduardo Feitosa
# 05 Ago 2025
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
        residues::Vector{Float64},
        Factors::Matrix{Float64},
        coords::Matrix{Int},
        sorted::Vector{Int},
        inlist::Matrix{Bool},
        listsize::Int,
        rbp_not_converged::Bool
    )
    
    @fastmath @inbounds for e in 1:num_reps

        # display("e = $e")

        # # display(sum(inlist))

        # display([residues coords'])

        # 1) Find largest residue  and coordenates
        if residues[1] == 0.0
            init_list_RBP!(Lq,Lr,Nc,nothing,nothing,newLr,Factors,inlist,
                                            residues,coords,listsize)
            if residues[1] == 0.0
                rbp_not_converged = false
                break
            end
        end
        cimax = coords[1,1]
        vjmax = coords[2,1]
        limax = coords[3,1]

        sorted = sort(coords[2,1:listsize])
        count_max = 1
        count = 1
        max_idx = 1
        for i=2:listsize
            if sorted[i] > 0
                if sorted[i] == sorted[i-1]
                    count += 1
                else
                    count = 1
                end
                if count > count_max
                    count_max = count
                    max_idx = i
                end
            end
        end
        # display(max_idx)
        if count_max > 1
            vjmax = sorted[max_idx]
        end

        # display(vjmax)

        Nvjmax = Nv[vjmax]
        for ci in Nvjmax
            li = LinearIndices(Lr)[ci,vjmax]
            # # display(li)
            # 2) Decay the RBP factor corresponding to the maximum residue
            Factors[li] *= decayfactor
            # 3) update check to node message Lr[cnmax,vnmax]
            Lr[li] = newLr[li]
            # 4) remove from list
            remove_list_VN(inlist,li,listsize,coords,residues)
        end 

        # # display([residues coords'])

        Ld = calc_Ld(vjmax,Nvjmax,Lf,Lr)
        bitvector[vjmax] = signbit(Ld)

        for ci in Nvjmax
            # 5) update Nv messages Lq[ci,vnmax]
            li = LinearIndices(Lq)[ci,vjmax]
            Lq[li] = tanh(0.5*(Ld - Lr[li]))
            # 6) calculate residues
            Nci = Nc[ci]
            for vj in Nci
                if vj ≠ vjmax
                    li = LinearIndices(Lr)[ci,vj]
                    newlr = calc_Lr(Nci,ci,vj,Lq)
                    li = LinearIndices(Lr)[ci,vj]
                    newLr[li] = newlr
                    residue = abs(newlr - Lr[li])*Factors[li]
                    remove_list_VN(inlist,li,listsize,coords,residues)
                    add_list_VN!(residue,residues,listsize,inlist,coords,li,ci,vj)
                end
            end
        end

        # # display([residues coords'])

    end

    return rbp_not_converged

end

function 
    remove_list_VN(
        inlist::Matrix{Bool},
        li::Int,
        listsize::Int,
        coords::Matrix{Int},
        residues::Vector{Float64}
    )

    if inlist[li]
        pos = 0
        for i = 1:listsize
            if coords[3,i] == li
                pos = i
                break
            end
        end
        if pos == 0
            throw(error("($(coords[1,i]),$(coords[2,i])) is registered as being on the list, but it's not."))
        end
        # remove from list
        inlist[li] = false
        # update the list
        for i in pos:listsize
            residues[i] = residues[i+1]
            coords[1,i] = coords[1,i+1]
            coords[2,i] = coords[2,i+1]
            coords[3,i] = coords[3,i+1]
        end
    end
end

function 
    add_list_VN!(
        residue::Float64,
        residues::Vector{Float64},
        listsize::Int,
        inlist::Matrix{Bool},
        coords::Matrix{Int},
        li::Int,
        ci::Int,
        vj::Int
    )

    if residue > residues[listsize]
        if residue ≥ residues[1]
            i = 1
        else
            d = listsize >> 1
            i = d
            while d > 1
                d >>= 1
                if residue ≥ residues[i]
                    i -= d
                else
                    i += d
                end
            end
            if residue < residues[i]
                i += 1
            end
        end

        last = coords[3,end-1]
        if last ≠ 0
            inlist[last] = false
        end
        inlist[li] = true

        for j=listsize:-1:i+1
            residues[j] = residues[j-1]
            coords[1,j] = coords[1,j-1]
            coords[2,j] = coords[2,j-1]
            coords[3,j] = coords[3,j-1]
        end
        coords[1,i] = ci
        coords[2,i] = vj
        coords[3,i] = li
        residues[i] = residue        
    end
end
