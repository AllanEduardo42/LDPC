################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Vertical update and MAP estimate of the LLR based Sum-Product Algorithm

### No Inf restriction ###
function 
    llr_vertical_update_and_MAP_crude!(
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

### Inf restriction ###

function 
    llr_vertical_update_and_MAP!(
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
        @inbounds Ld = ΔLf[node]
        for check in checks
            @inbounds @fastmath Ld += Lr[check,node]
        end
        if Ld < 0
            @inbounds d[node] = 1
        end
        for check in checks
            @inbounds @fastmath Lq[check,node] = Ld - Lr[check,node]
        end
    end
end


function 
    llr_init_q!(
        Lq::Matrix{<:AbstractFloat},
        ΔLf::Vector{<:AbstractFloat},
        nodes2checks::Vector{Vector{T}} where {T<:Integer}
    )

    node = 0
    for checks in nodes2checks
        node += 1
        for check in checks
            @inbounds Lq[check,node] = ΔLf[node]
        end
    end
end