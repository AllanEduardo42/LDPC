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
        Residues::Matrix{Float64},
        residues::Vector{Float64},
        Factors::Matrix{Float64},
        coords::Matrix{Int},
        inlist::Matrix{Bool},
        listsize::Int,
        listsize2::Int,
        rbp_not_converged::Bool
    )
    
    @fastmath @inbounds for e in 1:num_reps

        # display("e = $e")

        # display(sum(inlist))

        # display([residues coords'])

        # 1) Find largest residue and coordenates
        if residues[1] == 0.0
            for ci in eachindex(Nc)
                for vj in Nc[ci]
                    li = LinearIndices(Residues)[ci,vj]
                    residue = Residues[li]
                    add_residue!(inlist,residues,coords,residue,li,ci,vj,listsize,listsize2)
                end
            end
            if residues[1] == 0.0
                rbp_not_converged = false
                break
            end
        end
        cimax = coords[1,1]
        vjmax = coords[2,1]
        limax = coords[3,1]
              
        Nvjmax = Nv[vjmax]
        for ci in Nvjmax
            li = LinearIndices(Lr)[ci,vjmax]
            # 2) Decay the RBP factor corresponding to the maximum residue
            Factors[li] *= decayfactor
            # 3) update check to node message Lr[cnmax,vnmax]
            Lr[li] = newLr[li]
            # 4) remove from list
            Residues[li] = 0.0
            remove_list_VN(inlist,li,listsize,coords,residues)
        end

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
                    Residues[li] = residue
                    remove_list_VN(inlist,li,listsize,coords,residues)
                    add_residue!(inlist,residues,coords,residue,li,ci,vj,listsize,listsize2)
                end
            end
        end
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

    @inbounds if inlist[li]
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
