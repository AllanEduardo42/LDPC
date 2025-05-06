################################################################################
# Allan Eduardo Feitosa
# 31 Mar 2025
# RBP Sum-Product Algorithm with residual decaying factor

include("./RBP functions/findmaxnode.jl")
include("./RBP functions/calc_residue.jl")

function
    VN_RBP!(
        bitvector::Vector{Bool},
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        Nc::Vector{Vector{T}} where {T<:Integer},
        Nv::Vector{Vector{T}} where {T<:Integer},
        Lrj::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{<:AbstractFloat},Nothing},
        decayfactor::AbstractFloat,
        num_edges::Integer,
        Ms::Matrix{<:AbstractFloat},
        Factors::Vector{<:AbstractFloat},
        residues::Vector{<:AbstractFloat},
        bp_not_converged::Bool
    )

    @fastmath @inbounds for e in 1:num_edges

        # display("e = $e")

        vjmax = findmaxnode(residues)
        if vjmax == 0
            bp_not_converged = false
            break # i.e., BP has converged
        end
        Factors[vjmax] *= decayfactor
        residues[vjmax] = 0.0

        for ci in Nv[vjmax]
            li = LinearIndices(Ms)[ci,vjmax]
            Lr[li] = Ms[li]
        end        

        bitvector[vjmax] = update_Lq!(Lq,Lr,Lf,vjmax,Nv[vjmax],Lrj)

        for ci in Nv[vjmax]
            # calculate the new check to node messages
            Nci = Nc[ci]
            update_Lr!(Ms,Lq,ci,Nci,Lrj,signs,phi)
            for vj in Nci
                if vj ≠ vjmax
                    li = LinearIndices(Ms)[ci,vj]
                    residue = Ms[li] - Lr[li]
                    residues[vj] = abs(residue)*Factors[vj]
                end
            end
        end
    end

    return bp_not_converged

end

