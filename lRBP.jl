
################################################################################
# Allan Eduardo Feitosa
# 26 set 2024
# local RBP Sum-Product Algorithm using min-sum to calculate the residues

include("min_sum.jl")
include("check2node_llr.jl")
include("node2check_llr_and_MAP.jl")
include("RBP_vertical_update.jl")

function
    lRBP!(
        d::Vector{Bool},
        Lr::Matrix{<:AbstractFloat},
        max_coords::Vector{<:Integer},
        Lq::Matrix{<:AbstractFloat},
        ΔLf::Vector{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        nodes2checks::Vector{Vector{T}} where {T<:Integer},
        sn::Vector{Bool},
        Edges::Matrix{<:Integer},
        penalty::Matrix{<:AbstractFloat},
        penalty_factor::AbstractFloat,
        num_edges::Integer
    )

    e = 1
    while e <= num_edges

        (cmax,nmax) = max_coords
        penalty[cmax,nmax] *= penalty_factor
        Edges[cmax,nmax] += 1
        Lr[cmax,nmax] = check2node_llr_one_check_only!(
            view(Lq,cmax,:),
            checks2nodes[cmax],
            nmax,
            Lr[cmax,nmax]
        )
        max_residue = 0.0
        
        _checks = nodes2checks[nmax]
        for check in _checks
            if check ≠ cmax
                # vertical update of Lq[check,nmax], check ≠ cmax
                @inbounds Lq[check,nmax] =  RBP_vertical_update(
                                                _checks,
                                                check,
                                                ΔLf[nmax],
                                                view(Lr,:,nmax)
                                            )
                # find max ΔLr[node,check], node ≠ nmax, check ≠ cmax
                _nodes = checks2nodes[check]
                max_residue = min_sum_lRBP!(
                    max_coords,
                    max_residue,
                    penalty,
                    view(Lr,check,:),
                    view(Lq,check,:),
                    sn,
                    _nodes,
                    nmax,
                    check
                )
            end
        end

        if max_residue == 0.0
            break
        end

        e += 1

    end

    MAP!(
        d,
        nodes2checks,
        ΔLf,
        Lr
    )

end
