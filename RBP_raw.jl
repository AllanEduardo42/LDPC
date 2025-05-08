################################################################################
# Allan Eduardo Feitosa
# 06 Mai 2025
# RBP Algorithm with residual decaying factor (RAW)

include("./RBP functions/findmaxedge.jl")
include("./RBP functions/calc_residue.jl")

function
    RBP!(
        bitvector::Vector{Bool},
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        Nc::Vector{Vector{T}} where {T<:Integer},
        Nv::Vector{Vector{T}} where {T<:Integer},
        ::Nothing,
        ::Nothing,
        ::Nothing,
        decayfactor::AbstractFloat,
        num_edges::Integer,
        newLr::Matrix{<:AbstractFloat},
        Factors::Matrix{<:AbstractFloat},
        coords::Matrix{<:Integer},
        edges::Matrix{<:Integer},
        residues::Vector{<:AbstractFloat},
        ::Nothing,
        ::Nothing,
        ::Vector{<:Integer},
        relative::Bool,
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
            limax = coords[3,max_edge]
        end

        # 2) Decay the RBP factor corresponding to the maximum residue
        Factors[limax] *= decayfactor
        # 3) update check to node message Lr[cnmax,vnmax]
        Lr[limax] = newLr[limax]
        # 4) set maximum residue to zero
        residues[max_edge] = 0.0

        Nvjmax = Nv[vjmax]
        for ci in Nvjmax
            if ci ≠ cimax
                # 5) update Nv messages Lq[ci,vnmax]
                Lq[ci,vjmax] = calc_Lq(Nvjmax,ci,vjmax,Lr,Lf)
                # 6) calculate residues
                Nci = Nc[ci]
                for vj in Nci
                    if vj ≠ vjmax
                        li = LinearIndices(Lr)[ci,vj]
                        newlr = calc_Lr(Nci,ci,vj,Lq)
                        oldLr = Lr[li]
                        newLr[li] = newlr
                        residues[edges[li]] = calc_residue(newlr,oldLr,
                            Factors[li],relative,Lq[li],nothing)
                    end
                end
            end
        end
    end

    # 7) update bitvector
    for vj in eachindex(Nv)
        ci = Nv[vj][1]
        li = LinearIndices(Lr)[ci,vj]
        bitvector[vj] = signbit(Lr[li] + Lq[li])
    end

    return bp_not_converged
end