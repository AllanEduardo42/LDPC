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
        num_edges::Integer,
        ::Matrix{<:AbstractFloat},
        Factors::Vector{<:AbstractFloat},
        alpha::Vector{<:AbstractFloat},
        bp_not_converged::Bool
    )

    e = 0
    @fastmath @inbounds while e < num_edges

        # display("e = $e")

        # 1) Find largest residue  and coordenates
        ci = findmaxnode(alpha)
        # display(findmax(alpha))
        if ci == 0
            bp_not_converged = false
            break # i.e., RBP has converged
        end

        # 4) set maximum residue to zero
        alpha[ci] = 0.0

        # 2) Decay the RBP factor corresponding to the maximum residue
        # Factors[ci] *= decayfactor
        # 3) update check to node message Lr[cnmax,vnmax]
        Nci = Nc[ci]
        for vk in Nci
            e += 1
            Lr[ci,vk] = calc_Lr(Nci,ci,vk,Lq)
            Nvk = Nn[vk]
            for ca in Nvk
                if ca ≠ ci
                    # 5) update Nn messages Lq[ca,vnmax]
                    Lq[ca,vk] = calc_Lq(Nvk,ca,vk,Lr,Lf)   
                    # 6) calculate alpha
                    Nca = Nc[ca]
                    alpha[ca] = calc_alpha(Nca,ca,Lr,Lq,Factors)
                end
            end
        end
    end

    # 7) update bitvector
    for vb in eachindex(Nn)
        ca = Nn[vb][1]
        bitvector[vb] = signbit(Lr[ca,vb] + Lq[ca,vb])
    end

    return bp_not_converged
end

function calc_Lr(
    Nci::Vector{<:Integer},
    ci::Integer,
    vk::Integer,    
    Lq::Matrix{<:AbstractFloat}
)

    @fastmath @inbounds begin
        pLr = 1.0
        for vb in Nci
            if vb ≠ vk
                pLr *= tanh(0.5*Lq[ci,vb])
            end
        end
        return 2*atanh(pLr)
    end
end

function calc_Lq(
    Nvk::Vector{<:Integer},
    ci::Integer,
    vk::Integer,
    Lr::Matrix{<:AbstractFloat},
    Lf::Vector{<:AbstractFloat}
)

    @fastmath @inbounds begin
        Lq = Lf[vk]
        for ca in Nvk
            if ca ≠ ci
                Lq += Lr[ca,vk]
            end
        end
        return Lq
    end
end

function calc_residue(
    Nni::Vector{<:Integer},
    ni::Integer,
    nj::Integer,
    Lr::Matrix{<:AbstractFloat},
    Lq::Matrix{<:AbstractFloat},
    Factors::Vector{<:AbstractFloat},
)

    @inbounds begin
        residue = calc_Lr(Nni,ni,nj,Lq) - Lr[ni,nj]
        if isnan(residue)
            return 0.0
        else
            # @fastmath return abs(x)*Factors[ni]
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