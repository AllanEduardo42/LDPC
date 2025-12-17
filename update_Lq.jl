################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Vertical update and MAP estimate of the LLR based Sum-Product Algorithm (with
# "Inf" restriction)


include("tanhLq.jl")

############################ SPA USING LLRs METHOD #############################

function 
    calc_Ld(
        vj::Int,
        Nvj::Vector{Int},
        Lf::Vector{Float64},
        Lr::Matrix{Float64}
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
    Nvj::Vector{Int},
    ci::Int,
    vj::Int,
    Lr::Matrix{Float64},
    Lf::Vector{Float64}
)::Float64

    @fastmath @inbounds begin
        lq = Lf[vj]
        for ca in Nvj
            if ca ≠ ci
                lq += Lr[ca,vj]
            end
        end
        return tanh(0.5*lq)
    end
end

########################### SPA USING MKAY's METHOD ############################

function
    calc_Ld(
        r::Array{Float64,3},
        f::Matrix{Float64},
        ci::Int,
        vj::Int,
        Nvj::Vector{Int}
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