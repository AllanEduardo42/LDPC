################################################################################
# Allan Eduardo Feitosa
# 26 Feb 2025
# List-RBP Sum-Product Algorithm with residual decaying factor

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/calc_residues.jl")
include("./RBP functions/update_list.jl")
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
        listm1::Vector{<:Integer},
        listn1::Vector{<:Integer},
        listres2::Union{Vector{<:AbstractFloat},Nothing},
        listm2::Union{Vector{<:Integer},Nothing},
        listn2::Union{Vector{<:Integer},Nothing},
        inlist::Union{Matrix{<:Integer},Nothing}        
    )

    @fastmath @inbounds for e in 1:num_edges

        # display("e = $e")

        # if in test mode, store the values of the maximum residues
        if test
            all_max_res_alt[e] = listres1[1]     
        end
        
        # 1) get the largest residues coordenates if not clipped
        if listres1[1] == 0.0
            cnmax = rand(rng_rbp,1:length(cn2vn))
            vnmax = rand(rng_rbp,cn2vn[cnmax])
        else
            cnmax = listm1[1]
            vnmax = listn1[1]
        end

        # 2) Decay the RBP factor corresponding to the maximum residue
        lmax = decay!(cnmax,vnmax,Factors,decayfactor)

        # 3) update check to node message Lr[cnmax,vnmax]
        RBP_update_Lr!(lmax,Lr,Ms,cnmax,vnmax,cn2vn,Lq,Lrn,signs,phi)

        # 4) Remove max residue from the list and update the list
        remove_from_list!(lmax,listsizes[1],listres1,listm1,listn1,inlist,1) 

        # 5) update Ldn[vmax] and bitvector[vnmax]
        bitvector[vnmax] = update_Lq!(Lq,Lr,Lf[vnmax],vnmax,vn2cn,Lrn)

        # 6) update vn2cn messages Lq[vnmax,m], ∀m ≠ cnmax, and calculate residues
        calc_residues!(Factors,Ms,Lr,Lq,Lrn,signs,phi,cnmax,vnmax,cn2vn,
            vn2cn[vnmax],listres1,listm1,listn1,listres2,listm2,listn2,
            listsizes,inlist)
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
#             m = listm1[index]
#             n = listn1[index]
#             lmax = LinearIndices(Factors)[m,n]
#             remove_from_list!(lmax,listsizes[1],listres1,listm1,listn1,inlist,index)
#         end
#     end
# end

