################################################################################
# Allan Eduardo Feitosa
# 11 out 2024
# RBP Sum-Product Algorithm

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/RBP_set_zero_or_remove.jl")
include("./RBP functions/calc_residues.jl")
include("./RBP functions/findmaxcoords.jl")
include("./RBP functions/update_list.jl")
include("calc_syndrome.jl")

#RBP
function
    RBP!(
        address::Union{Matrix{<:Integer},Nothing},
        addressinv::Union{Matrix{<:Integer},Nothing},
        residues::Union{Vector{<:AbstractFloat},Nothing},
        bitvector::Vector{Bool},
        Lr::Matrix{<:AbstractFloat},
        Ms::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        Lrn::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},        
        phi::Union{Vector{<:AbstractFloat},Nothing},
        Factors::Matrix{<:AbstractFloat},
        decayfactor::AbstractFloat,
        num_edges::Integer,
        rng_rbp::AbstractRNG,
        listsize::Integer,
        listsize2::Integer,
        listres::Vector{<:AbstractFloat},
        listm::Vector{<:Integer},
        listn::Vector{<:Integer},
        listres2::Union{Vector{<:AbstractFloat},Nothing},
        listm2::Union{Vector{<:Integer},Nothing},
        listn2::Union{Vector{<:Integer},Nothing},
        inlist::Union{Matrix{<:Integer},Nothing},
        max_residues::Vector{<:AbstractFloat}
    )

    @inbounds @fastmath for e in 1:num_edges

        # display("e = $e")
        # display([listm listn listres])

        max_residues[e] = listres[1]       

        index = 1
        if listres[index] != 0
            cnmax = listm[index]
            vnmax = listn[index]
        else
            cnmax = rand(rng_rbp,1:length(cn2vn))
            vnmax = rand(rng_rbp,cn2vn[cnmax])
        end

        # display([cnmax vnmax])

        # 2) for optimization, get the linear indice of the maximum residue
        lmax = LinearIndices(Factors)[cnmax,vnmax]

        # 3) Decay the RBP factor corresponding to the maximum residue
        Factors[lmax] *= decayfactor

        # 4) update check to node message Lr[cnmax,vnmax]
        RBP_update_Lr!(lmax,cnmax,vnmax,cn2vn,Lq,Lr,Ms,Lrn,signs,phi)

        # 5) set maximum residue to zero or remove it from the list
        set_zero_or_remove!(addressinv,residues,lmax,listsize,listres,listm,
                            listn,inlist,index)

        # 6) update Ldn[vmax] and bitvector[vnmax]
        _, bitvector[vnmax] = update_Lq!(Lq,Lr,Lf[vnmax],vnmax,vn2cn,Lrn)

        # 7) update vn2cn messages Lq[vnmax,m], ∀m ≠ cnmax, and calculate residues
        leaf = true # suppose node vnmax is a leaf in the graph
        count_size = listsize2
        for m in vn2cn[vnmax]
            if m ≠ cnmax
                leaf = false # vnmax is not a leaf
                count_size = calc_residues!(addressinv,residues,Factors,Ms,Lr,Lq,Lrn,signs,
                               phi,vnmax,m,cn2vn,listres,listm,listn,listres2,
                               listm2,listn2,listsize,listsize2,count_size,inlist)
            end
        end

        # 8) if vnmax is a leaf in the graph
        if leaf
            count_size = calc_residues!(addressinv,residues,Factors,Ms,Lr,Lq,Lrn,signs,
                               phi,vnmax,cnmax,cn2vn,listres,listm,listn,listres2,
                               listm2,listn2,listsize,listsize2,count_size,inlist)
        end

        # display("list 1")
        # display([listm listn listres])
        # display("list 2")
        # display([listm2 listn2 listres2])

        # 9) Update list 1
        if listsize2 ≠ 0
            count = 0
            k = listsize
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
            for k in count_size:-1:1
                update_list!(inlist,listres,listm,listn,listres2[k],listm2[k],
                             listn2[k],listsize)
            end
            listres2 .*= 0.0
            listm2 .*= 0
            listn2 .*= 0
        end

        # display("list 1 updated")
        # display([listm listn listres])

        # 10) find maximum residue (only for the original RBP)

        findmaxcoords!(address,residues,listres,listm,listn)

    end
end

