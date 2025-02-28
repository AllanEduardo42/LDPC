################################################################################
# Allan Eduardo Feitosa
# 26 Feb 2025
# RBP Sum-Product Algorithm using the only local strategy with residual decaying
# factor

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/decay.jl")

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
        all_max_res_alt::Union{Vector{<:AbstractFloat},Nothing},
        test::Bool,
        rng_rbp::AbstractRNG,
        max_residue::Vector{<:AbstractFloat},
        max_coords::Vector{<:Integer}
    )

    @fastmath @inbounds for e in 1:num_edges

        # display("e = $e")

        # if in test mode, store the values of the maximum residues
        if test
            all_max_res_alt[e] = max_residue[1] 
        end 

        # 1) get the largest residues coordenates
        if max_residue[1] == 0.0
            cnmax = rand(rng_rbp,1:length(cn2vn))
            vnmax = rand(rng_rbp,cn2vn[cnmax])
        else
            cnmax = max_coords[1]
            vnmax = max_coords[2]
        end        

        # 2) Decay the RBP factor corresponding to the maximum residue
        lmax = decay!(cnmax,vnmax,Factors,decayfactor)

        # 3) update check to node message Lr[cnmax,vnmax]
        RBP_update_Lr!(lmax,Lr,Ms,cnmax,vnmax,cn2vn,Lq,Lrn,signs,phi)
        
        # 4) clear the max residue
        max_residue[1] = 0.0

        # 5) update Ldn[vmax] and bitvector[vnmax]
        bitvector[vnmax] = update_Lq!(Lq,Lr,Lf[vnmax],vnmax,vn2cn,Lrn)

        # 6) update vn2cn messages Lq[vnmax,m], ∀m ≠ cnmax, and calculate residues
        calc_residues!(max_residue,max_coords,Factors,Ms,Lr,Lq,Lrn,signs,phi,
            cnmax,vnmax,cn2vn,vn2cn[vnmax])
    end
end

