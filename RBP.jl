
include("MAP.jl")
include("RBP_vertical_update.jl")
include("RBP_findmax_ΔLr.jl")

function
    RBP!(
        d::Vector{Bool},
        Lr::Matrix{<:AbstractFloat},
        max_coords::Vector{<:Integer},
        max_residue::AbstractFloat,
        Lq::Matrix{<:AbstractFloat},
        ΔLf::Vector{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        nodes2checks::Vector{Vector{T}} where {T<:Integer},
        Lrn::Vector{<:AbstractFloat},
    )

    for m in eachindex(checks2nodes)

        if max_residue == 0.0
            break
        end

        (cmax,nmax) = max_coords

        # update the message with largest residue
        @inbounds @fastmath Lr[cmax,nmax] += max_residue
        max_residue = 0.0
        
        checks2nodes_nmax = nodes2checks[nmax]
        for check in checks2nodes_nmax
            if check ≠ cmax
                # vertical update of Lq[check,nmax], check ≠ cmax
                @inbounds Lq[check,nmax] =  RBP_vertical_update(
                                                checks2nodes_nmax,
                                                check,
                                                ΔLf[nmax],
                                                view(Lr,:,nmax)
                                            )
                # find max ΔLr[node,check], node ≠ nmax, check ≠ cmax
                _checks2nodes = checks2nodes[check]
                max_residue =   RBP_findmax_ΔLr!(
                                    _checks2nodes,
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
    end

    MAP!(
        d,
        nodes2checks,
        ΔLf,
        Lr
    )

end