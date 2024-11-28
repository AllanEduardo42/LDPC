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
    cn2vn::Vector{Vector{T}} where {T<:Integer},
    alpha::AbstractFloat,
    alpha2::AbstractFloat
    )

    s = false
    minL = INFFLOAT
    minL2 = INFFLOAT
    max_idx = 0
    @fastmath @inbounds for n in cn2vn[m]
        β, signs[n], s = abs_sign!(Lq[n,m],s)
        if β < minL
            max_idx = n
            minL, minL2 = β, minL
        elseif β < minL2
            minL2 = β
        end
    end

    @fastmath @inbounds for n in cn2vn[m]
        if n == max_idx #(pick the second least Lq)
            Ms[m,n] = (alpha - alpha2*(signs[n] ⊻ s))*minL2
        else
            Ms[m,n] = (alpha - alpha2*(signs[n] ⊻ s))*minL
        end
    end

end

