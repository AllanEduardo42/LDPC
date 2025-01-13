################################################################################
# Allan Eduardo Feitosa
# 13 jan 2025
# RBP Sum-Product Algorithm using the only local strategy

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/RBP_set_zero_or_remove.jl")
include("./RBP functions/calc_residues.jl")
include("./RBP functions/find_local_maxresidue.jl")
include("./RBP functions/update_list.jl")

#Local-RBP
function
    local_RBP!(
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
        decay::AbstractFloat,
        num_edges::Integer,
        Ldn::Vector{<:AbstractFloat},
        rng_sample::AbstractRNG,
    )

    @fastmath @inbounds for e in 1:num_edges

        # 1) if maximum residue is zero, the RBP has converged
        if maxresidue == 0
            break
        end
        
        # 2) get the check and node of the maximum residue
        cnmax = maxcoords[1]
        vnmax = maxcoords[2]

        # 3) verify if the list was not updated
        if cnmax == 0 || maxresidue == -1
            cnmax = rand(rng_sample,1:length(cn2vn))
            vnmax = rand(rng_sample,cn2vn[cnmax])
        end

        # 4) for optimization, get the linear indice of the maximum residue
        lmax = LinearIndices(Factors)[cnmax,vnmax]

        # 5) Decay the RBP factor corresponding to the maximum residue
        Factors[lmax] *= decay

        # 6) update check to node message Lr[cnmax,vnmax]
        RBP_update_Lr!(lmax,cnmax,vnmax,cn2vn,Lq,Lr,Ms,Lrn,signs,phi)

        # 7) set maximum residue to zero or remove it from the list
        maxresidue = 0

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
                maxresidue = find_local_maxresidue!(Factors,Ms,Lr,Lq,Lrn,signs,
                               phi,vnmax,m,cn2vn,maxcoords)
            end
        end

        # 10) if vnmax is a leaf in the graph, triggers the random selection of 
        #     a check in 3)
        if leaf
            maxresidue = -1
        end

    end
end

