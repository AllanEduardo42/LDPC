
################################################################################
# Allan Eduardo Feitosa
# 23 set 2024
# RBP Sum-Product Algorithm using min-sum to calculate the residues

include("min_sum.jl")
include("llr_horizontal_update.jl")
include("MAP.jl")
include("RBP_vertical_update.jl")

function
    RBP_R!(
        d::Vector{Bool},
        Lr::Matrix{<:AbstractFloat},
        max_coords::Vector{<:Integer},
        Lq::Matrix{<:AbstractFloat},
        ΔLf::Vector{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        nodes2checks::Vector{Vector{T}} where {T<:Integer},
        sn::Vector{Int8},
        R::Matrix{<:AbstractFloat},
        Edges::Matrix{<:Integer}
    )

    for m in 1:EDGES

        max_residue = find_max_residue_coords!(max_coords,R,checks2nodes)

        if max_residue == 0.0 # if RBP has converged
            break
        end

        (cmax,nmax) = max_coords
        Edges[cmax,nmax] += 1
        Lr[cmax,nmax] = llr_horizontal_update_one_check_only!(
            view(Lq,cmax,:),
            checks2nodes[cmax],
            nmax,
            Lr[cmax,nmax]
        )
        # we don't use the check-to-node message correspoding to (cmax,nmax) anymore.
        # Thus we make the largest residue equal to zero:
        R[cmax,nmax] = 0.0
        
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
                # if any new residue estimate is larger than the previously estimated maximum 
                # residue than update the value of max_residue and max_coords.
                min_sum_RBP_R!(
                    R,
                    view(Lr,check,:),
                    view(Lq,check,:),
                    sn,
                    _nodes,
                    nmax,
                    check
                )
            end
        end

    end

    MAP!(
        d,
        nodes2checks,
        ΔLf,
        Lr
    )

end

function
    find_max_residue_coords!(
        max_coords::Vector{<:Integer},
        R::Matrix{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer}
    )

    check = 0
    max_residue = 0
    x = 0.0
    for nodes in checks2nodes
        check += 1
        for node in nodes
            x = R[check,node]
            if x > max_residue
                max_residue = x
                max_coords[1] = check
                max_coords[2] = node
            end
        end
    end

    return max_residue
end