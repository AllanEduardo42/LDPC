################################################################################
# Allan Eduardo Feitosa
# 22 Feb 2024
# RBP Sum-Product Algorithm with residual decaying factor

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/calc_residues.jl")
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
        address::Union{Matrix{<:Integer},Nothing},
        addressinv::Union{Matrix{<:Integer},Nothing},
        residues::Union{Vector{<:AbstractFloat},Nothing}
    )

    @inbounds @fastmath for e in 1:num_edges

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

        # 4) set maximum residue to zero or remove it from the list
        residues[max_edge] = 0.0

        # 5) update Ldn[vmax] and bitvector[vnmax]
        bitvector[vnmax] = update_Lq!(Lq,Lr,Lf[vnmax],vnmax,vn2cn,Lrn)
        # Obs: this code update Lq[vnmax,cnmax]: no difference in performance 
        #      was detected. The original implementation, which doesn't update
        #      Lq[vnmax,cnmax], is as follows (steps 1, 2 and 3):
        # step 1) Ldn[vnmax], nl = calc_Ld(vnmax,vn2cn,Lf[vnmax],Lr)
        # step 2) bitvector[vnmax] = signbit(Ldn[vnmax]
        # step 3) see below.

        # 6) update vn2cn messages Lq[vnmax,m], ∀m ≠ cnmax, and calculate residues
        checks = vn2cn[vnmax]
        if length(checks) > 1
            for m in checks
                if m ≠ cnmax
                    # step 3) Lq[vnmax,m] = Ldn[vnmax] - Lr[nl+m]        
                    calc_residues!(addressinv,residues,Factors,Ms,Lr,Lq,Lrn,signs,
                                phi,vnmax,m,cn2vn)
                end
            end
        else # if vnmax is a leaf in the graph
            calc_residues!(addressinv,residues,Factors,Ms,Lr,Lq,Lrn,signs,
                           phi,vnmax,m,cn2vn)
        end
    end
end

