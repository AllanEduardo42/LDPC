################################################################################
# Allan Eduardo Feitosa
# 22 Feb 2024
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
        Ldn::Vector{<:AbstractFloat},
        Ms::Matrix{<:AbstractFloat},
        Factors::Matrix{<:AbstractFloat},        
        all_max_res_alt::Union{Vector{<:AbstractFloat},Nothing},
        test::Bool,
        rng_rbp::AbstractRNG,
        listsize1::Integer,
        listsize2::Integer,
        listres::Vector{<:AbstractFloat},
        listm::Vector{<:Integer},
        listn::Vector{<:Integer},
        listres2::Union{Vector{<:AbstractFloat},Nothing},
        listm2::Union{Vector{<:Integer},Nothing},
        listn2::Union{Vector{<:Integer},Nothing},
        inlist::Union{Matrix{<:Integer},Nothing}        
    )

    @inbounds @fastmath for e in 1:num_edges

        # display("e = $e")

        if test
            all_max_res_alt[e] = listres[1]     
        end
        
        # 1) get the largest residues coordenates if not clipped
        new_listsize2 = listsize2
        index = 1
        if listsize1 > 1
            search_index = true
            i = 0
            while search_index
                i += 1
                if abs(listres[index]) < CLIP
                    search_index = false
                    index = i
                else
                    new_listsize2 += 1
                    m = listm[index]
                    n = listn[index]
                    lmax = LinearIndices(Factors)[m,n]
                    remove_from_list!(lmax,listsize1,listres,listm,listn,inlist,index)
                end
            end
        end
        if listres[index] != 0
            cnmax = listm[index]
            vnmax = listn[index]
        else
            cnmax = rand(rng_rbp,1:length(cn2vn))
            vnmax = rand(rng_rbp,cn2vn[cnmax])
        end        

        # 2) Decay the RBP factor corresponding to the maximum residue
        lmax = decay!(cnmax,vnmax,Factors,decayfactor)

        # 3) update check to node message Lr[cnmax,vnmax]
        RBP_update_Lr!(lmax,Lr,Ms,cnmax,vnmax,cn2vn,Lq,Lrn,signs,phi)

        # 4) Remove max residue from the list and update the list
        remove_from_list!(lmax,listsize1,listres,listm,listn,inlist,index) 

        # 5) update Ldn[vmax] and bitvector[vnmax]
        Ldn[vnmax], nl = calc_Ld(vnmax,vn2cn,Lf[vnmax],Lr)
        bitvector[vnmax] = signbit(Ldn[vnmax])

        # 6) update vn2cn messages Lq[vnmax,m], ∀m ≠ cnmax, and calculate residues
        leaf = true # suppose node vnmax is a leaf in the graph
        for m in vn2cn[vnmax]
            if m ≠ cnmax
                leaf = false # vnmax is not a leaf
                Lq[vnmax,m] = Ldn[vnmax] - Lr[nl+m]
                new_listsize2 = calc_residues!(Factors,Ms,Lr,Lq,Lrn,signs,
                               phi,vnmax,m,cn2vn,listres,listm,listn,listres2,
                               listm2,listn2,listsize1,listsize2,new_listsize2,inlist)
            end
        end

        # 7) if vnmax is a leaf in the graph
        if leaf
            Lq[vnmax,cnmax] = Ldn[vnmax] - Lr[nl+cnmax]
            new_listsize2 = calc_residues!(Factors,Ms,Lr,Lq,Lrn,signs,
                               phi,vnmax,cnmax,cn2vn,listres,listm,listn,listres2,
                               listm2,listn2,listsize1,listsize2,new_listsize2,inlist)
        end

        # 8) update list
        if listsize2 ≠ 0
            count = 0
            k = listsize1
            while count < listsize2-1
                count += 1
                k -= 1
                m = listm[k]
                if m != 0
                    listres[k] = 0.0
                    n = listn[k]
                    inlist[m,n] = false
                end
            end
            for k in new_listsize2:-1:1
                update_list!(inlist,listres,listm,listn,listres2[k],listm2[k],
                             listn2[k],listsize1)
            end
            listres2 .*= 0.0
            listm2 .*= 0
            listn2 .*= 0
        end     

    end
end

