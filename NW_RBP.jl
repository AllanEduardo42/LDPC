################################################################################
# Allan Eduardo Feitosa
# 07 Apr 2025
# RBP and List-RBP Algorithm with residual decaying factor

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/findmaxedge.jl")
include("./RBP functions/update_lists.jl")
include("./RBP functions/remove_residue!.jl")

function
    NW_RBP!(
        bitvector::Vector{Bool},
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        Nc::Vector{Vector{T}} where {T<:Integer},
        Nn::Vector{Vector{T}} where {T<:Integer},
        ::Union{Vector{<:AbstractFloat},Nothing},
        ::Union{Vector{Bool},Nothing},
        ::Union{Vector{<:AbstractFloat},Nothing},
        decayfactor::AbstractFloat,
        M::Integer,
        ::Matrix{<:AbstractFloat},
        Factors::Vector{<:AbstractFloat},
        alpha::Vector{<:AbstractFloat},
        bp_not_converged::Bool
    )

    @fastmath @inbounds for m = 1:M

        # display("e = $e")

        # 1) Find largest residue  and coordenates
        cimax = findmaxnode(alpha)
        # display(findmax(alpha))
        if cimax == 0
            bp_not_converged = false
            break # i.e., RBP has converged
        end

        # 4) set maximum residue to zero
        alpha[cimax] = 0.0

        # 2) Decay the RBP factor corresponding to the maximum residue
        # Factors[cimax] *= decayfactor
        # 3) update check to node message Lr[cnmax,vnmax]
        Ncimax = Nc[cimax]
        for vj in Ncimax
            Lr[cimax,vj] = calc_Lr(Ncimax,cimax,vj,Lq)
            Nvj = Nn[vj]
            for ci in Nvj
                if ci ≠ cimax
                    # 5) update Nn messages Lq[ci,vnmax]
                    Lq[ci,vj] = calc_Lq(Nvj,ci,vj,Lr,Lf)   
                    # 6) calculate alpha
                    Nci = Nc[ci]
                    alpha[ci] = calc_alpha(Nci,ci,Lr,Lq,Factors)
                end
            end
        end
    end

    # 7) update bitvector
    for vb in eachindex(Nn)
        ci = Nn[vb][1]
        bitvector[vb] = signbit(Lr[ci,vb] + Lq[ci,vb])
    end

    return bp_not_converged
end

function calc_residue(
    Nci::Vector{<:Integer},
    ci::Integer,
    vj::Integer,
    Lr::Matrix{<:AbstractFloat},
    Lq::Matrix{<:AbstractFloat},
    Factors::Vector{<:AbstractFloat},
)

    @inbounds begin
        residue = calc_Lr(Nci,ci,vj,Lq) - Lr[ci,vj]
        if isnan(residue)
            return 0.0
        else
            # @fastmath return abs(x)*Factors[ci]
            @fastmath return abs(residue)
        end
    end
end

function calc_alpha(
    Nci::Vector{<:Integer},
    ci::Integer,
    Lr::Matrix{<:AbstractFloat},
    Lq::Matrix{<:AbstractFloat},
    Factors::Vector{<:AbstractFloat}
)

    alpha = 0.0
    @fastmath @inbounds for vb in Nci
        residue = calc_residue(Nci,ci,vb,Lr,Lq,Factors)
        if residue > alpha
            alpha = residue
        end
    end

    return alpha
end