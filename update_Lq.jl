################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Vertical update and MAP estimate of the LLR based Sum-Product Algorithm (with
# "Inf" restriction)

############################ SPA USING LLRs METHOD #############################

function 
    calc_Ld(
        vj::Integer,
        Nvj::Vector{<:Integer},
        Lf::Vector{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat}
    )

    @fastmath @inbounds begin
        Ld = Lf[vj]
        for ci in Nvj
            Ld += Lr[ci,vj]
        end
    end
    
    return Ld

end

######################### SPA USING LLRs METHOD NO OPT #########################

function calc_Lq(
    Nvj::Vector{<:Integer},
    ci::Integer,
    vj::Integer,
    Lr::Matrix{<:AbstractFloat},
    Lf::Vector{<:AbstractFloat}
)

    @fastmath @inbounds begin
        Lq = Lf[vj]
        for ca in Nvj
            if ca ≠ ci
                Lq += Lr[ca,vj]
            end
        end
        return Lq
    end
end

########################### SPA USING MKAY's METHOD ############################

function
    calc_Ld(
        r::Array{<:AbstractFloat,3},
        f::Matrix{<:AbstractFloat},
        ci::Integer,
        vj::Integer,
        Nvj::Vector{<:Integer}
    )

    @fastmath @inbounds begin
        Ld1 = f[vj,1]
        Ld2 = f[vj,2]
        for ca in Nvj
            if ca ≠ ci
                Ld1 *= r[ca,vj,1]
                Ld2 *= r[ca,vj,2]
            end
        end
        return Ld1, Ld2
    end
end