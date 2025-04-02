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
    list_RBP_genius!(
        bitvector::Vector{Bool},
        cword::Vector{Bool},
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
        all_max_res_alt::Union{Vector{<:AbstractFloat},Nothing},
        test::Bool,
        rng_rbp::AbstractRNG,
        listsizes::Vector{<:Integer},
        listres1::Vector{<:AbstractFloat},
        indices_res1::Vector{<:Integer},
        listres2::Union{Vector{<:AbstractFloat},Nothing},
        indices_res2::Union{Vector{<:Integer},Nothing},
        inlist::Union{Matrix{<:Integer},Nothing}        
    )

    listsize1m1 = listsizes[1]-1
    difflistsizes = listsizes[1] - listsizes[2]+1

    @fastmath @inbounds for e in 1:num_edges

        # display("e = $e")
        # display([listres1 indices_res1])

        # if in test mode, store the values of the maximum residues
        if test
            all_max_res_alt[e] = listres1[1]     
        end

        wrong = true
        index = 1
        lmax = 1
        cnmax = 1
        vnmax = 1
        newbitvector = false
        while wrong && index ≤ listsizes[1]
            # 1) get the largest residues coordenates if not clipped
            if listres1[1] == 0.0
                cnmax = rand(rng_rbp,1:length(cn2vn))
                vnmax = rand(rng_rbp,cn2vn[cnmax])
                lmax = LinearIndices(Factors)[cnmax,vnmax]
            else
                lmax = indices_res1[1]
                ci = CartesianIndices(Factors)[lmax]
                cnmax, vnmax = ci[1], ci[2]
            end
            Lr_old = Lr[lmax]
            Lr[lmax] = Ms[lmax]
            newbitvector = update_Lq!(Lq,Lr,Lf[vnmax],vnmax,vn2cn,Lrn)

            if newbitvector == cword[vnmax]
                wrong = false
            else
                Lr[lmax] = Lr_old
                index += 1
            end
        end

        # 2) Decay the RBP factor corresponding to the maximum residue
        Factors[lmax] *= decayfactor

        # 4) Remove max residue from the list and update the list
        remove_from_list!(lmax,listsizes[1],listres1,indices_res1,inlist,1)

        # 5) update vn2cn messages Lq[vnmax,m] and bitvector[vnmax]
        bitvector[vnmax] = newbitvector

        # 6) calculate local residues
        new_listsize2 = calc_local_residues_list!(Lq,Lr,cn2vn,vn2cn,Lrn,signs,phi,Ms,
            Factors,listsizes,listres1,indices_res1,listres2,indices_res2,inlist,
            listsizes[2],cnmax,vnmax)
        # update list 1 
        update_list1!(listres1,indices_res1,listres2,indices_res2,listsizes,
            new_listsize2,inlist,listsize1m1,difflistsizes)

        # clear list 2
        listres2 .*= 0.0
        indices_res2 .*= 0    
    end
end