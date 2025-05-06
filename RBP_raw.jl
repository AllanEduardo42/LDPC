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
        Nv::Vector{Vector{T}} where {T<:Integer},
        ::Union{Vector{<:AbstractFloat},Nothing},
        ::Union{Vector{Bool},Nothing},
        ::Union{Vector{<:AbstractFloat},Nothing},
        decayfactor::AbstractFloat,
        num_edges::Integer,
        ::Matrix{<:AbstractFloat},
        Factors::Matrix{<:AbstractFloat},
        coords::Matrix{<:Integer},
        edges::Matrix{<:Integer},
        residues::Vector{<:AbstractFloat},
        ::Nothing,
        ::Nothing,
        ::Vector{<:Integer},
        ::Bool,
        bp_not_converged::Bool
    )

    @fastmath @inbounds for e in 1:num_edges

        # display("e = $e")

        # 1) Find largest residue  and coordenates
        max_edge, _ = findmaxedge(residues,nothing)
        if max_edge == 0
            bp_not_converged = false
            break # i.e., RBP has converged
        else
            cimax = coords[1,max_edge]
            vjmax = coords[2,max_edge]
        end

        # 2) Decay the RBP factor corresponding to the maximum residue
        Factors[cimax,vjmax] *= decayfactor
        # 3) update check to node message Lr[cnmax,vnmax]
        Ncimax = Nc[cimax]
        Lr[cimax,vjmax] = calc_Lr(Ncimax,cimax,vjmax,Lq)
        # 4) set maximum residue to zero
        residues[max_edge] = 0.0

        Nvjmax = Nv[vjmax]
        for ci in Nvjmax
            if ci ≠ cimax
                # 5) update Nv messages Lq[ci,vnmax]
                Lq[ci,vjmax] = calc_Lq(Nvjmax,ci,vjmax,Lr,Lf)
                Nci = Nc[ci]    
                # 6) calculate residues
                for vj in Nci
                    if vj ≠ vjmax
                        residues[edges[ci,vj]] = calc_residue(Nci,ci,vj,Lr,Lq,Factors)
                    end
                end
            end
        end
    end

    # 7) update bitvector
    for vj in eachindex(Nv)
        ci = Nv[vj][1]
        bitvector[vj] = signbit(Lr[ci,vj] + Lq[ci,vj])
    end

    return bp_not_converged
end

function calc_residue(
    Nci::Vector{<:Integer},
    ci::Integer,
    vj::Integer,
    Lr::Matrix{<:AbstractFloat},
    Lq::Matrix{<:AbstractFloat},
    Factors::Matrix{<:AbstractFloat},
)

    @inbounds begin
        residue = calc_Lr(Nci,ci,vj,Lq) - Lr[ci,vj] # fastmath deactivated
        if isnan(residue)
            return 0.0
        else
            @fastmath return abs(residue)*Factors[ci,vj]
        end
    end
end