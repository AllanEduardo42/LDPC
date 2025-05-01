################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Vertical update and MAP estimate of the LLR based Sum-Product Algorithm (with
# "Inf" restriction)

############################ SPA USING LLRs METHOD #############################

function
    update_Lq!(
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lf::AbstractFloat,
        n::Integer,
        cns::Vector{<:Integer},
        ::Vector{<:AbstractFloat}
    )

    Ld = calc_Ld(n,cns,Lf,Lr)
    @fastmath @inbounds for m in cns
        li = LinearIndices(Lq)[m,n]
        Lq[li] = Ld - Lr[li]
    end

    return signbit(Ld)
end

function 
    calc_Ld(
        n::Integer,
        cns::Vector{<:Integer},
        Ld::AbstractFloat,
        Lr::Matrix{<:AbstractFloat}
    )

    @fastmath @inbounds for m in cns
        Ld += Lr[m,n]
    end
    
    return Ld

end

######################### SPA USING LLRs METHOD NO OPT #########################

function
    update_Lq!(
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lf::AbstractFloat,
        n::Integer,
        cns::Vector{<:Integer},
        ::Nothing
    )
    
    m = 0
    @inbounds for outer m in cns
        Lq[m,n] = Lf[n]
        for m2 in cns
            if m2 ≠ m
                Lq[m,n] += Lr[m2,n]
            end
        end
    end

    @inbounds Ld = Lq[m,n] + Lr[m,n]
    return signbit(Ld)
end

########################### SPA USING MKAY's METHOD ############################

function
    update_Lq!(
        q::Array{<:AbstractFloat,3},
        r::Array{<:AbstractFloat,3},
        f::Vector{<:AbstractFloat},
        n::Integer,
        cns::Vector{<:Integer}
    )

    @inbounds for m in cns
        Ld1 = f[1]
        Ld2 = f[2]
        for m2 in cns
            if m2 ≠ m
                Ld1 *= r[m2,n,1]
                Ld2 *= r[m2,n,2]
            end
        end
        a = Ld1 + Ld2
        q[m,n,1] = Ld1/a
        q[m,n,2] = Ld2/a
    end
end