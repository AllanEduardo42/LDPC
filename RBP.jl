################################################################################
# Allan Eduardo Feitosa
# 11 out 2024
# RBP Sum-Product Algorithm

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/RBP_set_zero_or_remove.jl")
include("./RBP functions/calc_residues.jl")
include("./RBP functions/findmaxcoords.jl")
include("./RBP functions/update_list.jl")

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
        Ldn::Vector{<:AbstractFloat},
        rng_sample::AbstractRNG,
        listsize::Integer,
        listsize2::Integer,
        listres::Vector{<:AbstractFloat},
        listm::Vector{<:Integer},
        listn::Vector{<:Integer},
        listres2::Union{Vector{<:AbstractFloat},Nothing},
        listm2::Union{Vector{<:Integer},Nothing},
        listn2::Union{Vector{<:Integer},Nothing},
        inlist::Union{Matrix{<:Integer},Nothing}
    )

    @fastmath @inbounds for e in 1:num_edges

        # 1) if maximum residue is zero, the RBP has converged
        # if listres[1] == 0
        #     break
        # end
        
        # 2) get the check and node of the maximum residue
        # cnmax = listm[1]
        # vnmax = listn[1]

        # 3) verify if the list was not updated
        if listres[1] != 0
            cnmax = listm[1]
            vnmax = listn[1]
        else
            cnmax = rand(rng_sample,1:length(cn2vn))
            vnmax = rand(rng_sample,cn2vn[cnmax])
        end

        # 4) for optimization, get the linear indice of the maximum residue
        lmax = LinearIndices(Factors)[cnmax,vnmax]

        # 5) Decay the RBP factor corresponding to the maximum residue
        Factors[lmax] *= decayfactor

        # 6) update check to node message Lr[cnmax,vnmax]
        RBP_update_Lr!(lmax,cnmax,vnmax,cn2vn,Lq,Lr,Ms,Lrn,signs,phi)

        # 7) set maximum residue to zero or remove it from the list
        set_zero_or_remove!(addressinv,residues,lmax,listsize,listres,listm,
                            listn,inlist)

        # 8) update Ldn[vmax] and bitvector[vnmax]
        Ldn[vnmax] = Lf[vnmax]
        nl = LinearIndices(Lr)[1,vnmax]-1
        for m in vn2cn[vnmax]
            Ldn[vnmax] += Lr[nl+m]
            bitvector[vnmax] = signbit(Ldn[vnmax])
        end

        # 9) update vn2cn messages Lq[vnmax,m], ∀m ≠ cnmax, and calculate residues
        leaf = true # suppose node vnmax is a leaf in the graph
        for m in vn2cn[vnmax]
            if m ≠ cnmax
                leaf = false # vnmax is not a leaf
                Lq[vnmax,m] = Ldn[vnmax] - Lr[nl+m]
                calc_residues!(addressinv,residues,Factors,Ms,Lr,Lq,Lrn,signs,
                               phi,vnmax,m,cn2vn,listres,listm,listn,listres2,
                               listm2,listn2,listsize,listsize2,inlist)
            end
        end

        # 10) if vnmax is a leaf in the graph, triggers the random selection of 
        #     a check in 3)
        if leaf && listsize == 1
            listres[1] = -1
        end

        # 11) Update list 2
        if listsize2 ≠ 0
            for k in listsize2:-1:1
                update_list!(inlist,listres,listm,listn,listres2[k],listm2[k],
                             listn2[k],listsize)
            end
            listres2 .*= 0.0
            listm2 .*= 0
            listn2 .*= 0
        end

        # 12) find maximum residue (only for the original RBP)

        findmaxcoords!(address,residues,listres,listm,listn)

    end
end

