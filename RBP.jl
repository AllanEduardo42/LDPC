
include("MAP.jl")
include("RBP_vertical_update.jl")
include("RBP_findmax_ΔLr.jl")

function
    RBP!(
        d::Vector{Bool},
        Lr::Matrix{<:AbstractFloat},
        max_coords::Vector{<:Integer},
        Lq::Matrix{<:AbstractFloat},
        ΔLf::Vector{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        nodes2checks::Vector{Vector{T}} where {T<:Integer},
        Lrn::Vector{<:AbstractFloat},
        R::Matrix{<:AbstractFloat}
    )

    max_residue = 0.0
    check = 0
    for nodes in checks2nodes
        check+=1
        max_residue =   RBP_findmax_ΔLr!(
                                    nodes,
                                    Lrn,
                                    view(Lq,check,:),
                                    view(Lr,check,:),
                                    check,
                                    0,
                                    max_residue,
                                    max_coords
                                )
        # calc_residues!(
        #     nodes,
        #     Lrn,
        #     view(Lq,check,:),
        #     view(Lr,check,:),
        #     check,
        #     0,
        #     R
        # )
    end

    # max_residue = find_max_residue(R,checks2nodes,max_coords)

    for m in 1:1531

        (cmax,nmax) = max_coords

        # update the message with largest residue
        @inbounds @fastmath Lr[cmax,nmax] += max_residue
        # R[cmax,nmax] = 0.0
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
                # calc_residues!(
                #     _nodes,
                #     Lrn,
                #     view(Lq,check,:),
                #     view(Lr,check,:),
                #     check,
                #     nmax,
                #     R
                #     )
                max_residue =   RBP_findmax_ΔLr!(
                                    _nodes,
                                    Lrn,
                                    view(Lq,check,:),
                                    view(Lr,check,:),
                                    check,
                                    nmax,
                                    max_residue,
                                    max_coords
                                )
            end
        end
        # max_residue = find_max_residue(R,checks2nodes,max_coords)

        if max_residue == 0.0
            break
        end
    end

    MAP!(
        d,
        nodes2checks,
        ΔLf,
        Lr
    )

end

function find_max_residue(R,checks2nodes,max_coords)

    check = 0
    max_residue = 0
    for nodes in checks2nodes
        check += 1
        for node in nodes
            if abs(R[check,node]) > abs(max_residue)
                max_residue = R[check,node]
                max_coords[1] = check
                max_coords[2] = node
            end
        end
    end

    return max_residue
end