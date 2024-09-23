################################################################################
# Allan Eduardo Feitosa
# 23 set 2024
# Horizontal and vertical update of the LLR based Sum-Product Algorithm

# These functions must be used if one wants to have no control on divergent values of
# of Lr and Lq (i.e., Lr[.,.] = Inf or Lq[.,.] = Inf).
# These functions avoids producing NaN values.

function 
    llr_horizontal_update_Inf!(
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer}
    )

    check = 0
    for nodes in checks2nodes
        check += 1
        for node in nodes
            @inbounds Lr[check,node] = 1.0
            for n in nodes
                if n ≠ node
                    @inbounds @fastmath Lr[check,node] *= tanh(0.5*Lq[check,n])
                end
            end
            @inbounds @fastmath Lr[check,node] = 2*atanh(Lr[check,node])
        end
    end
end

function 
    llr_vertical_update_and_MAP_Inf!(
        Lq::Matrix{<:AbstractFloat},
        d::Vector{Bool},
        Lr::Matrix{<:AbstractFloat},
        ΔLf::Vector{<:AbstractFloat},
        nodes2checks::Vector{Vector{T}} where {T<:Integer}
    )

    d .*= false
    node = 0
    for checks in nodes2checks
        node += 1
        for check in checks
            @inbounds Lq[check,node] = ΔLf[node]
            for c in checks
                if c ≠ check
                    @inbounds @fastmath Lq[check,node] += Lr[c,node]
                end
            end
        end
        @inbounds Ld = ΔLf[node]
        for check in checks
            Ld += Lr[check,node]
        end
        if Ld < 0
            @inbounds d[node] = 1
        end
    end
end