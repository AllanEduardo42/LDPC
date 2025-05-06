################################################################################
# Allan Eduardo Feitosa
# 07 Apr 2025
# RBP and List-RBP Algorithm with residual decaying factor

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/findmaxedge.jl")
include("./RBP functions/update_lists.jl")
include("./RBP functions/remove_residue!.jl")

function
    RBP_raw!(
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
        Factors::Matrix{<:AbstractFloat},
        coords::Matrix{<:Integer},
        indices::Matrix{<:Integer},
        residues::Vector{<:AbstractFloat},
        local_residues::Union{Vector{<:AbstractFloat},Nothing},
        ::Union{Matrix{<:Integer},Nothing},
        ::Vector{<:Integer},
        ::Bool,
        bp_not_converged::Bool
    )

    @fastmath @inbounds for e in 1:num_edges

        # display("e = $e")

        # 1) Find largest residue  and coordenates
        max_edge, _ = findmaxedge(residues,local_residues)
        if max_edge == 0
            bp_not_converged = false
            break # i.e., RBP has converged
        else
            ci = coords[1,max_edge]
            vj = coords[2,max_edge]
        end

        # 2) Decay the RBP factor corresponding to the maximum residue
        # Factors[limax] *= decayfactor
        # 3) update check to node message Lr[cnmax,vnmax]
        Nci = Nc[ci]
        Lr[ci,vj] = calc_Lr(Nci,ci,vj,Lq)
        # 4) set maximum residue to zero
        residues[max_edge] = 0.0

        Nvj = Nn[vj]
        for ca in Nvj
            if ca ≠ ci
                # 5) update Nn messages Lq[ca,vnmax]
                Lq[ca,vj] = calc_Lq(Nvj,ca,vj,Lr,Lf)
                Nca = Nc[ca]    
                # 6) calculate residues
                for vb in Nca
                    if vb ≠ vj
                        residues[indices[ca,vb]] = calc_residue(Nca,ca,vb,Lr,Lq,Factors)
                    end
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
    vj::Integer,    
    Lq::Matrix{<:AbstractFloat}
)

    @fastmath @inbounds begin
        pLr = 1.0
        for vb in Nci
            if vb ≠ vj
                pLr *= tanh(0.5*Lq[ci,vb])
            end
        end
        return 2*atanh(pLr)
    end
end

function calc_Lq(
    Nvj::Vector{<:Integer},
    ci::Integer,
    vj::Integer,
    Lr::Matrix{<:AbstractFloat},
    Lf::Vector{<:AbstractFloat}
)

    @fastmath @inbounds begin
        Lq = Lf[vj]
        for ca in Nvj
            if ca ≠ ci
                Lq += Lr[ca,vj]
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
    Factors::Matrix{<:AbstractFloat},
)

    @inbounds begin
        residue = calc_Lr(Nni,ni,nj,Lq) - Lr[ni,nj]
        if isnan(residue)
            return 0.0
        else
            # @fastmath return abs(x)*Factors[ca,vb]
            @fastmath return abs(residue)
        end
    end
end