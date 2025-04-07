################################################################################
# Allan Eduardo Feitosa
# 3 Mar 2025
# RBP Sum-Product Algorithm using the only local strategy with residual decaying
# factor

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/calc_local_residues.jl")

function
    local_RBP!(
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
        max_residue::Vector{<:AbstractFloat},
        maxcoords::Vector{<:Integer}
    )

    @fastmath @inbounds for e in 1:num_edges

        # display("e = $e")
        # display(max_residue)

        # 1) get the largest residues coordenates
        if max_residue[1] == 0.0
            cnmax = rand(rng_rbp,1:length(cn2vn))
            vnmax = rand(rng_rbp,cn2vn[cnmax])
        else
            cnmax = maxcoords[1]
            vnmax = maxcoords[2]
        end        

        # 2) Decay the RBP factor corresponding to the maximum residue
        Factors[cnmax,vnmax] *= decayfactor

        # 3) update check to node message Lr[cnmax,vnmax]
        # RBP_update_Lr!(lmax,Lr,Ms,cnmax,vnmax,cn2vn,Lq,Lrn,signs,phi)
        Lr[cnmax,vnmax] = Ms[cnmax,vnmax]
        
        # 4) clear the max residue
        max_residue[1] = 0.0
        # max_residue[2] = 0.0

        # 5) update vn2cn messages Lq[vnmax,m] and bitvector[vnmax]
        bitvector[vnmax] = update_Lq!(Lq,Lr,Lf[vnmax],vnmax,vn2cn[vnmax],Lrn)

        # 6) calculate residues
       calc_local_residues_local!(Lq,Lr,cn2vn,vn2cn[vnmax],Lrn,signs,phi,Ms,Factors,
            max_residue,maxcoords,cnmax,vnmax)
    
        # update list
        if max_residue[1] < max_residue[3]
            max_residue[1], max_residue[3] = max_residue[3], max_residue[1] 
            maxcoords[1], maxcoords[5] = maxcoords[5], maxcoords[1] 
            maxcoords[2], maxcoords[6] = maxcoords[6], maxcoords[2]      
        else
            max_residue[3] = max_residue[2]
            maxcoords[5] = maxcoords[3]
            maxcoords[6] = maxcoords[4]
        end
    end
end

