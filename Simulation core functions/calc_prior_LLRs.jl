################################################################################
# Allan Eduardo Feitosa
# 19 set 2024
# Function to calculate the prior log-likelihood ratios

function calc_prior_LLRs!(
    prior_LLR::Vector{Float64},
    twoZc::Int,
    signal::Vector{Float64},
    variance::Float64,
    rayleigh::Bool,
    fading::Union{Vector{Float64},Nothing},
    bptype::String
)

    @fastmath @inbounds begin
        if rayleigh
            for i in eachindex(signal)
                prior_LLR[twoZc+i] = -2*signal[i]*fading[i]/variance
            end
        else
            for i in eachindex(signal)
                prior_LLR[twoZc+i] = -2*signal[i]/variance
            end
        end
        if bptype == "TABL"
            for i in eachindex(signal)
                prior_LLR[twoZc+i] *= SIZE_PER_RANGE
            end
        end
    end

end

# MKAY
function calc_prior_LLRs!(
    f::Matrix{Float64},
    twoZc::Int,
    signal::Vector{Float64},
    variance::Float64
)

    @fastmath @inbounds begin

        k = 1/sqrt(2*π*variance)

        for i in eachindex(signal)
            f[twoZc+i,1] = k*exp(-(signal[i]+1)^2/(2*variance))
            f[twoZc+i,2] = k*exp(-(signal[i]-1)^2/(2*variance))
        end

        normalize!(f)
    end

end

function normalize!(f::AbstractMatrix{Float64})
    
    @inbounds @fastmath for n in axes(f,1)
        α = f[n,1] + f[n,2]
        f[n,1] = f[n,1]/α
        f[n,2] = f[n,2]/α
    end

end