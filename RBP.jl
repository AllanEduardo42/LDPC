################################################################################
# Allan Eduardo Feitosa
# 26 Feb 2025
# RBP Sum-Product Algorithm with residual decaying factor

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/calc_local_residues.jl")
include("./RBP functions/findmaxedge.jl")
include("./RBP functions/decay.jl")

function
    RBP!(
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
        address::Matrix{<:Integer},
        addressinv::Matrix{<:Integer},
        residues::Vector{<:AbstractFloat},
        thres::AbstractFloat
    )

    @fastmath @inbounds for e in 1:num_edges

        # display("e = $e")

        # 1) Find largest residue  and coordenates
        max_edge, maxresidue = findmaxedge(residues)
        if max_edge == 0
            break # i.e., RBP has converged
        else
            cnmax = address[1,max_edge]
            vnmax = address[2,max_edge]
        end
        # if in test mode, store the values of the maximum residues
        if test
            all_max_res_alt[e] = maxresidue
        end

        # 2) Decay the RBP factor corresponding to the maximum residue
        lmax = decay!(cnmax,vnmax,Factors,decayfactor)

        # 3) update check to node message Lr[cnmax,vnmax]
        RBP_update_Lr!(lmax,Lr,Ms,cnmax,vnmax,cn2vn,Lq,Lrn,signs,phi)

        # 4) set maximum residue to zero
        residues[max_edge] = 0.0

        # 5) update vn2cn messages Lq[vnmax,m] and bitvector[vnmax]
        bitvector[vnmax] = update_Lq!(Lq,Lr,Lf[vnmax],vnmax,vn2cn,Lrn)

        # 6) calculate residues
        calc_local_residues!(Lq,Lr,cn2vn,vn2cn,Lrn,signs,phi,Ms,Factors,addressinv,
            residues,cnmax,vnmax,thres)
    end
end

