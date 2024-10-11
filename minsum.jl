################################################################################
# Allan Eduardo Feitosa
# 22 set 2024
# Horizontal update of the LMs based MIN SUM Algorithm

function abs_sign!(Lq::AbstractFloat,s::Bool)
    sig = signbit(Lq)
    @fastmath ab = abs(Lq)
    return ab, sig, s ⊻ sig
end

function minsum!(                           
    Lq::Matrix{<:AbstractFloat},
    Ms::Matrix{<:AbstractFloat},
    signs::Vector{Bool},
    m::Integer,
    cn2vn::Vector{Vector{T}} where {T<:Integer}   
    )

    s = false
    minL = INF
    minL2 = INF
    max_idx = 0
    for n in cn2vn[m]
        @inbounds β, signs[n], s = abs_sign!(Lq[n,m],s)
        if @fastmath β < minL
            max_idx = n
            minL, minL2 = β, minL
        elseif @fastmath β < minL2
            minL2 = β
        end
    end

    for n in cn2vn[m]
        if n == max_idx #(pick the second least Lq)
            @fastmath @inbounds Ms[m,n] = (1 - 2*(signs[n] ⊻ s))*minL2
        else
            @fastmath @inbounds Ms[m,n] = (1 - 2*(signs[n] ⊻ s))*minL
        end
    end

end

