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
        largest_res::Vector{<:AbstractFloat},
        largestcoords::Vector{<:Integer},
        largestcoords_alt::Vector{<:Integer},
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
        rng_sample::AbstractRNG
    )

    largest_res_alt = largest_res[2]

    @fastmath @inbounds for e in 1:num_edges

        if largest_res[1] == 0
            cnmax = rand(rng_sample,1:length(cn2vn))
            vnmax = rand(rng_sample,cn2vn[cnmax])
        else
            largest_res[1] = 0
            cnmax = largestcoords[1]
            vnmax = largestcoords[2]
        end

        # 4) for optimization, get the linear indice of the maximum residue
        lmax = LinearIndices(Factors)[cnmax,vnmax]

        # 5) Decay the RBP factor corresponding to the maximum residue
        Factors[lmax] *= decay

        # 6) update check to node message Lr[cnmax,vnmax]
        RBP_update_Lr!(lmax,cnmax,vnmax,cn2vn,Lq,Lr,Ms,Lrn,signs,phi)        

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
                find_local_maxresidue!(largest_res,Factors,Ms,Lr,Lq,
                    Lrn,signs,phi,vnmax,m,cn2vn,largestcoords)
            end
        end

        # 10) if vnmax is a leaf in the graph, triggers the random selection of 
        #     a check in 3)
        if leaf
            Lq[vnmax,cnmax] = Ldn[vnmax] - Lr[nl+cnmax]
            find_local_maxresidue!(largest_res,Factors,Ms,Lr,Lq,
                Lrn,signs,phi,vnmax,cnmax,cn2vn,largestcoords)
        end
        if largest_res[1] < largest_res_alt
            largestcoords[1], largestcoords_alt[1] = largestcoords_alt[1], largestcoords[1]
            largestcoords[2], largestcoords_alt[2] = largestcoords_alt[2], largestcoords[2]
            largest_res[1], largest_res_alt = largest_res_alt, largest_res[1]            
        else
            largest_res_alt = largest_res[2]
            largestcoords_alt[1] = largestcoords[3]
            largestcoords_alt[2] = largestcoords[4]
        end
    end
end

