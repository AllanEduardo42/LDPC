################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Vertical update and MAP estimate of the LLR based Sum-Product Algorithm (with
# "Inf" restriction)

########################### SPA USING MKAY's METHOD ############################

function
    update_Lq!(
        q::Array{<:AbstractFloat,3},
        r::Array{<:AbstractFloat,3},
        f::Vector{<:AbstractFloat},
        n::Integer,
        vn2cn::Vector{Vector{T}} where {T<:Integer}
    )

    Ld = zeros(2)
    for m in vn2cn[n]
        Ld .= f
        for m2 in vn2cn[n]
            if m2 ≠ m
                @fastmath @inbounds Ld[1] *= r[m2,n,1]
                @fastmath @inbounds Ld[2] *= r[m2,n,2]
            end
        end
        @fastmath a = sum(Ld)
        @fastmath @inbounds q[n,m,1] = Ld[1]/a
        @fastmath @inbounds q[n,m,2] = Ld[2]/a
        @fastmath @inbounds Ld[1] *= r[m,n,1]
        @fastmath @inbounds Ld[2] *= r[m,n,2]
    end

    return @fastmath @inbounds signbit(Ld[1]-Ld[2])

end

############################ SPA USING LLRs METHOD #############################

function
    update_Lq!(
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lf::AbstractFloat,
        n::Integer,
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        ::Vector{<:AbstractFloat}
    )

    Ld = calc_Ld(n,vn2cn,Lf,Lr)
    for m in vn2cn[n]
        @fastmath @inbounds Lq[n,m] = Ld - Lr[m,n]
    end

    return Ld, signbit(Ld)
end

function 
    calc_Ld(
        n::Integer,
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        Lf::AbstractFloat,
        Lr::Matrix{<:AbstractFloat}
    )
    Ld = Lf
    for m in vn2cn[n]
        @fastmath @inbounds Ld += Lr[m,n]
    end
    
    return Ld

end

######################### SPA USING LLRs METHOD VER 2 ##########################

function
    update_Lq!(
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lf::AbstractFloat,
        n::Integer,
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        ::Nothing
    )

    for m in vn2cn[n]
        @inbounds Lq[n,m] = Lf[n]
        for m2 in vn2cn[n]
            if m2 ≠ m
                @fastmath @inbounds Lq[n,m] += Lr[m2,n]
            end
        end
    end
    @inbounds Ld = Lf[n]
    for m in vn2cn[n]
        @fastmath @inbounds Ld += Lr[m,n]
    end

    return Ld, signbit(Ld)
end