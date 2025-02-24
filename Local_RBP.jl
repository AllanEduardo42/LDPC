################################################################################
# Allan Eduardo Feitosa
# 13 jan 2025
# RBP Sum-Product Algorithm using the only local strategy with residual decaying
# factor

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/find_local_maxresidue.jl")
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
        max_coords::Vector{<:Integer},
        max_coords_alt::Vector{<:Integer}
    )

    max_residue_alt = max_residue[2]

    @fastmath @inbounds for e in 1:num_edges

        # display("e = $e")

        # if in test mode, store the values of the maximum residues
        if test
            all_max_res_alt[e] = max_residue[1] 
        end 

        # 1) get the largest residues coordenates
        if max_residue[1] == 0
            cnmax = rand(rng_rbp,1:length(cn2vn))
            vnmax = rand(rng_rbp,cn2vn[cnmax])
        else
            max_residue[1] = 0
            cnmax = max_coords[1]
            vnmax = max_coords[2]
        end        

        # 2) Decay the RBP factor corresponding to the maximum residue
        lmax = decay!(cnmax,vnmax,Factors,decayfactor)

        # 3) update check to node message Lr[cnmax,vnmax]
        RBP_update_Lr!(lmax,Lr,Ms,cnmax,vnmax,cn2vn,Lq,Lrn,signs,phi)    

        # 4) update Ldn[vmax] and bitvector[vnmax]
        bitvector[vnmax] = update_Lq!(Lq,Lr,Lf[vnmax],vnmax,vn2cn,Lrn)
        # Obs: this code update Lq[vnmax,cnmax]: no difference in performance 
        #      was detected. The original implementation, which doesn't update
        #      Lq[vnmax,cnmax], is as follows (steps 1, 2 and 3):
        # step 1) Ldn[vnmax], nl = calc_Ld(vnmax,vn2cn,Lf[vnmax],Lr)
        # step 2) bitvector[vnmax] = signbit(Ldn[vnmax]
        # step 3) see below.

        # 5) update vn2cn messages Lq[vnmax,m], ∀m ≠ cnmax, and calculate residues
        checks = vn2cn[vnmax]
        if length(checks) > 1
            for m in checks
                if m ≠ cnmax
                    # step 3) Lq[vnmax,m] = Ldn[vnmax] - Lr[nl+m]        
                    find_local_maxresidue!(max_residue,Factors,Ms,Lr,Lq,
                        Lrn,signs,phi,vnmax,m,cn2vn,max_coords)
                end
            end
        else # if vnmax is a leaf in the graph
            find_local_maxresidue!(max_residue,Factors,Ms,Lr,Lq,Lrn,signs,phi,
                vnmax,m,cn2vn,max_coords)
        end

        # 6) update list
        if max_residue[1] < max_residue_alt
            max_coords[1], max_coords_alt[1] = max_coords_alt[1], max_coords[1]
            max_coords[2], max_coords_alt[2] = max_coords_alt[2], max_coords[2]
            max_residue[1], max_residue_alt = max_residue_alt, max_residue[1]            
        else
            max_residue_alt = max_residue[2]
            max_coords_alt[1] = max_coords[3]
            max_coords_alt[2] = max_coords[4]
        end
    end
end

