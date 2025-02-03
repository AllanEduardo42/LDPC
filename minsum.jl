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
    minL = INFFLOAT
    minL2 = INFFLOAT
    max_idx = cn2vn[m][1]
    ml = LinearIndices(Lq)[1,m]-1
    @fastmath @inbounds for n in cn2vn[m]
        β, signs[n], s = abs_sign!(Lq[ml+n],s)
        if β < minL
            max_idx = n
            minL, minL2 = β, minL
        elseif β < minL2
            minL2 = β
        end
    end
    minL *= ALPHA
    minL2 *= ALPHA

    @fastmath @inbounds for n in cn2vn[m]
        if signs[n] ⊻ s
            Ms[m,n] = -minL
        else
            Ms[m,n] = minL
        end
    end
    @fastmath @inbounds if signbit(Ms[m,max_idx])
        Ms[m,max_idx] = -minL2
    else
        Ms[m,max_idx] = minL2
    end

end

