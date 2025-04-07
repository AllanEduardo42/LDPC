################################################################################
# Allan Eduardo Feitosa
# 31 Mar 2025
# RBP Sum-Product Algorithm with residual decaying factor

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/calc_residue.jl")
include("./RBP functions/findmaxnode.jl")

function
    NRBP!(
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
        Factors::Vector{<:AbstractFloat},
        residues::Vector{<:AbstractFloat},
        relative::Bool,
        rbp_not_converged::Bool
    )

    @fastmath @inbounds for e = 1:num_edges

        # display("e = $e")

        # display(sort(residues,rev=true))

        vnmax = findmaxnode(residues)
        if vnmax == 0
            rbp_not_converged = false
            break # i.e., RBP has converged
        end

        residues[vnmax] = 0.0

        Factors[vnmax] *= decayfactor

        for m in vn2cn[vnmax]
            Lr[m,vnmax] = Ms[m,vnmax]
        end        

        bitvector[vnmax] = update_Lq!(Lq,Lr,Lf[vnmax],vnmax,vn2cn,Lrn)

        for m in vn2cn[vnmax] 
            # calculate the new check to node messages
            update_Lr!(Ms,Lq,m,cn2vn,Lrn,signs,phi)
            for n in cn2vn[m]
                if n ≠ vnmax
                    li = LinearIndices(Ms)[m,n]
                    residues[n] = calc_residue(Ms,Lr,Factors,Lrn,Lq,li,relative)   
                end
            end
        end
    end

    return rbp_not_converged
end

