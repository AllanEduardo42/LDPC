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
                @inbounds @fastmath Ld[1] *= r[m2,n,1]
                @inbounds @fastmath Ld[2] *= r[m2,n,2]
            end
        end
        @fastmath a = sum(Ld)
        @inbounds @fastmath q[n,m,1] = Ld[1]/a
        @inbounds @fastmath q[n,m,2] = Ld[2]/a
        @inbounds @fastmath Ld[1] *= r[m,n,1]
        @inbounds @fastmath Ld[2] *= r[m,n,2]
    end

    return signbit(Ld[1]-Ld[2])

end

############################ SPA USING LLRs METHOD #############################

function
    update_Lq!(
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lf::AbstractFloat,
        n::Integer,
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        Lrn::Vector{<:AbstractFloat}
    )

    Ld = calc_Ld(n,vn2cn,Lf,Lr)
    for m in vn2cn[n]
        @inbounds @fastmath Lq[n,m] = Ld - Lr[m,n]
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
        @inbounds @fastmath Ld += Lr[m,n]
    end
    
    return Ld

end

function 
    MAP!(
        d::Vector{Bool},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        Lf::Vector{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat}
    )
    
    Ld = 0.0
    for n in eachindex(vn2cn)
        Ld = calc_Ld(n,vn2cn,Lf[n],Lr)
        @inbounds d[n] = signbit(Ld)
    end

end

######################### SPA USING LLRs METHOD VER 2 ##########################

function
    update_Lq!(
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lf::AbstractFloat,
        n::Integer,
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        Lrn::Nothing
    )

    for m in vn2cn[n]
        @inbounds Lq[n,m] = Lf[n]
        for m2 in vn2cn[n]
            if m2 ≠ m
                @inbounds @fastmath Lq[n,m] += Lr[m2,n]
            end
        end
    end
    @inbounds Ld = Lf[n]
    for m in vn2cn[n]
        Ld += Lr[m,n]
    end

    return Ld, signbit(Ld)
end