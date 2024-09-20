################################################################################
# Allan Eduardo Feitosa
# 19 set 2024
# Vertical update of the RBP Sum-Product Algorithm

function
    RBP_vertical_update(
        nodes2checks_nmax::Vector{<:Integer},
        check::Integer,
        Lq::AbstractFloat,
        Lr::AbstractVector{<:AbstractFloat}
    )
    
    for c in nodes2checks_nmax
        if c != check
            @inbounds @fastmath Lq += Lr[c]
        end
    end

    return Lq
end