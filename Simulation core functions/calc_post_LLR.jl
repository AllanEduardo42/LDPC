################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Calculate the posterior Log Likelihood Ratio (LLR) for SPA

############################ SPA USING LLRs METHOD #############################

function calc_post_LLR(              
    vj::Int,
    Nvj::Vector{Int},
    Lf::Vector{Float64},
    C2V::Matrix{Float64}
)

    # begin
    @fastmath @inbounds begin
        Ld = Lf[vj]
        for ci in Nvj
            Ld += C2V[ci,vj]
        end
    end
    
    return Ld

end

########################### SPA USING MKAY's METHOD ############################

function calc_post_LLR(
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