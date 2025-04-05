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
    vns::Vector{<:Integer}
    )

    @fastmath @inbounds begin

        s = false
        minL = INFFLOAT
        minL2 = INFFLOAT
        max_idx = vns[1]
        for n in vns
            β, signs[n], s = abs_sign!(Lq[m,n],s)
            if β < minL
                max_idx = n
                minL, minL2 = β, minL
            elseif β < minL2
                minL2 = β
            end
        end
        minL *= ALPHA
        minL2 *= ALPHA

        for n in vns
            if signs[n] ⊻ s
                Ms[m,n] = -minL
            else
                Ms[m,n] = minL
            end
        end
        if signbit(Ms[m,max_idx])
            Ms[m,max_idx] = -minL2
        else
            Ms[m,max_idx] = minL2
        end

    end

end

