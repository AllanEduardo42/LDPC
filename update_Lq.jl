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
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        ::Vector{<:AbstractFloat}
    )

    Ld, nl = calc_Ld(n,vn2cn,Lf,Lr)
    # for m in vn2cn[n]
    @fastmath @inbounds for m in vn2cn[n]
        Lq[n,m] = Ld - Lr[nl+m]
    end

    return signbit(Ld)
end

function 
    calc_Ld(
        n::Integer,
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        Lf::AbstractFloat,
        Lr::Matrix{<:AbstractFloat}
    )

    # for m in vn2cn[n]
    nl = LinearIndices(Lr)[1,n]-1
    @fastmath @inbounds for m in vn2cn[n]
        Lf += Lr[nl+m]
    end
    
    return Lf, nl

end

######################### SPA USING LLRs METHOD NO OPT #########################

function
    update_Lq!(
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lf::AbstractFloat,
        n::Integer,
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        ::Nothing
    )
    
    m = 0
    @inbounds for outer m in vn2cn[n]
        Lq[n,m] = Lf[n]
        for m2 in vn2cn[n]
            if m2 ≠ m
                Lq[n,m] += Lr[m2,n]
            end
        end
    end

    @inbounds Ld = Lq[n,m] + Lr[m,n]
    return signbit(Ld)
end

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
    @inbounds for m in vn2cn[n]
        Ld .= f
        @inbounds for m2 in vn2cn[n]
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