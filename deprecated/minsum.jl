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
    Lr::Matrix{<:AbstractFloat},
    signs::Vector{Bool},
    ci::Integer,
    Nci::Vector{<:Integer}
    )

    @fastmath @inbounds begin

        s = false
        minL = INFFLOAT
        minL2 = INFFLOAT
        max_idx = Nci[1]
        for vj in Nci
            β, signs[vj], s = abs_sign!(Lq[ci,vj],s)
            if β < minL
                max_idx = vj
                minL, minL2 = β, minL
            elseif β < minL2
                minL2 = β
            end
        end
        minL *= ALPHA
        minL2 *= ALPHA

        for vj in Nci
            if signs[vj] ⊻ s
                Lr[ci,vj] = -minL
            else
                Lr[ci,vj] = minL
            end
        end
        if signbit(Lr[ci,max_idx])
            Lr[ci,max_idx] = -minL2
        else
            Lr[ci,max_idx] = minL2
        end

    end

end

