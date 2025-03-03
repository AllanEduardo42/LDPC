################################################################################
# Allan Eduardo Feitosa
# 26 Feb 2025
# List-RBP Sum-Product Algorithm with residual decaying factor

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/calc_residue.jl")
include("./RBP functions/update_lists.jl")
include("./RBP functions/decay.jl")
include("./RBP functions/remove_from_list.jl")

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

        # 2) Decay the RBP factor corresponding to the maximum residue
        # lmax = decay!(cnmax,vnmax,Factors,decayfactor)
        Factors[lmax] *= decayfactor

        # 3) update check to node message Lr[cnmax,vnmax]
        RBP_update_Lr!(lmax,Lr,Ms,cnmax,vnmax,cn2vn,Lq,Lrn,signs,phi)

        # 4) Remove max residue from the list and update the list
        remove_from_list!(lmax,listsizes[1],listres1,indices_res1,inlist,1)

        # 5) update vn2cn messages Lq[vnmax,m] and bitvector[vnmax]
        bitvector[vnmax] = update_Lq!(Lq,Lr,Lf[vnmax],vnmax,vn2cn,Lrn)

        # 6) calculate residues
        new_listsize2 = listsizes[2]
        for m in vn2cn[vnmax]
            if m ≠ cnmax    
                # calculate the new check to node messages
                update_Lr!(Ms,Lq,m,cn2vn,Lrn,signs,phi)
                # calculate the residues
                for n in cn2vn[m]
                    if n ≠ vnmax
                        residue, li = calc_residue(Ms,Lr,Factors,Lrn,Lq,m,n)
                        new_listsize2 = update_list2!(listres1,indices_res1,
                                            listres2,indices_res2,listsizes,
                                            new_listsize2,inlist,li,residue)
                    end
                end
            end
        end
        # update list 1 
        update_list1!(listres1,indices_res1,listres2,indices_res2,listsizes,
            new_listsize2,inlist,listsize1m1,difflistsizes)

        # clear list 2
        listres2 .*= 0.0
        indices_res2 .*= 0    
    end
end

# if listsizes[1] > 1
#     search_index = true
#     i = 0
#     while search_index
#         i += 1
#         if abs(listres1[index]) < CLIP
#             search_index = false
#             index = i
#         else
#             if listsizes[2] == 1
#                 new_listsize2 += 1
#             end
#             m = indices_res1[index]
#             n = listn1[index]
#             lmax = LinearIndices(Factors)[m,n]
#             remove_from_list!(lmax,listsizes[1],listres1,indices_res1,listn1,inlist,index)
#         end
#     end
# end

