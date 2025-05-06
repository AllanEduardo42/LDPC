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
        Lf::Vector{<:AbstractFloat},
        vj::Integer,
        Nvj::Vector{<:Integer},
        ::Vector{<:AbstractFloat}
    )

    Ld = calc_Ld(vj,Nvj,Lf,Lr)
    @fastmath @inbounds for ci in Nvj
        li = LinearIndices(Lq)[ci,vj]
        Lq[li] = Ld - Lr[li]
    end

    return signbit(Ld)
end

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

function
    update_Lq!(
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        vj::Integer,
        Nvj::Vector{<:Integer},
        ::Nothing
    )
    
    ci = 0
    @fastmath @inbounds begin 
        for outer ci in Nvj
            Lq[ci,vj] = calc_Lq(Nvj,ci,vj,Lr,Lf)
        end
        # get the last ci of the loop iteration to calc Ld
        Ld = Lq[ci,vj] + Lr[ci,vj]
    end
    return signbit(Ld)
end

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
    update_Lq!(
        q::Array{<:AbstractFloat,3},
        r::Array{<:AbstractFloat,3},
        f::Vector{<:AbstractFloat},
        vj::Integer,
        Nvj::Vector{<:Integer}
    )

    @inbounds for ci in Nvj
        Ld1 = f[1]
        Ld2 = f[2]
        for ca in Nvj
            if ca ≠ ci
                Ld1 *= r[ca,vj,1]
                Ld2 *= r[ca,vj,2]
            end
        end
        a = Ld1 + Ld2
        q[ci,vj,1] = Ld1/a
        q[ci,vj,2] = Ld2/a
    end
end