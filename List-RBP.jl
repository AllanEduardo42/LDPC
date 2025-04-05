################################################################################
# Allan Eduardo Feitosa
# 26 Feb 2025
# List-RBP Sum-Product Algorithm with residual decaying factor

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/calc_residue.jl")
include("./RBP functions/update_lists.jl")
include("./RBP functions/decay.jl")
include("./RBP functions/remove_from_list.jl")
include("./RBP functions//calc_local_residues.jl")

function
    list_RBP!(
        bitvector::Vector{Bool},
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        Lrn::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{<:AbstractFloat},Nothing},
        decayfactor::AbstractFloat,
        num_edges::Integer,
        Ms::Matrix{<:AbstractFloat},
        Factors::Matrix{<:AbstractFloat},
        rng_rbp::AbstractRNG,
        listsizes::Vector{<:Integer},
        listres1::Vector{<:AbstractFloat},
        coords1::Matrix{<:Integer},
        listres2::Union{Vector{<:AbstractFloat},Nothing},
        coords2::Union{Matrix{<:Integer},Nothing},
        inlist::Union{Matrix{<:Integer},Nothing}        
    )

    listsize1m1 = listsizes[1]-1
    difflistsizes = listsizes[1] - listsizes[2]+1

    @fastmath @inbounds for e in 1:num_edges

        # display("e = $e")
        
        # 1) get the largest residues coordenates if not clipped
        if listres1[1] == 0.0
            cnmax = rand(rng_rbp,1:length(cn2vn))
            vnmax = rand(rng_rbp,cn2vn[cnmax])
        else
            cnmax = coords1[1,1]
            vnmax = coords1[2,1]
        end

        # 2) Decay the RBP factor corresponding to the maximum residue
        # lmax = decay!(cnmax,vnmax,Factors,decayfactor)
        Factors[cnmax,vnmax] *= decayfactor

        # 3) update check to node message Lr[cnmax,vnmax]
        # RBP_update_Lr!(lmax,Lr,Ms,cnmax,vnmax,cn2vn,Lq,Lrn,signs,phi)
        Lr[cnmax,vnmax] = Ms[cnmax,vnmax]

        # 4) Remove max residue from the list and update the list
        remove_from_list!(cnmax,vnmax,listsizes[1],listres1,coords1,inlist,1)

        # 5) update vn2cn messages Lq[vnmax,m] and bitvector[vnmax]
        bitvector[vnmax] = update_Lq!(Lq,Lr,Lf[vnmax],vnmax,vn2cn[vnmax],Lrn)

        # 6) calculate local residues
        new_listsize2 = calc_local_residues_list!(Lq,Lr,cn2vn,vn2cn[vnmax],Lrn,
            signs,phi,Ms,Factors,listsizes,listres1,coords1,listres2,coords2,
            inlist,listsizes[2],cnmax,vnmax)
        # update list 1 
        update_list1!(listres1,coords1,listres2,coords2,listsizes,new_listsize2,
            inlist,listsize1m1,difflistsizes)

        # clear list 2
        listres2 .*= 0.0
        coords2 .*= 0    
    end
end
