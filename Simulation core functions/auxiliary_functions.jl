################################################################################
# Allan Eduardo Feitosa
# 26 mai 2025
# Auxiliary functions used in sincore.jl

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

function init_V2C!(
    V2C::Matrix{Float64},
    prior_LLRs::Vector{Float64},
    Nv::Vector{Vector{Int}},
    msum_factor::Union{Float64,Nothing}
)
    
    @inbounds for vj in eachindex(Nv)
        for ci in Nv[vj]
            V2C[ci,vj] = tanh_V2C(prior_LLRs[vj],0.0,msum_factor)
        end
    end
end

function received_signal!(
    signal::Vector{Float64},
    cword::Vector{Bool},
    G::Int,
    twoLs::Int,
    σ::Float64,
    rgn::AbstractRNG,
    rayleigh::Bool,
    fading::Union{Vector{Float64},Nothing},
    x1::Union{Vector{Float64},Nothing},
    x2::Union{Vector{Float64},Nothing}
)

    @inbounds @fastmath begin
        randn!(rgn,signal)          # put the noise in the vector 'signal'
        lmul!(σ,signal)             # multiply by the standard deviation
        if rayleigh
            randn!(rgn,x1)
            randn!(rgn,x2)
            for g in 1:G
                u = 2*cword[twoLs+g] - 1
                fading[g] = sqrt(0.5*(x1[g]^2 + x2[g]^2))  
                signal[g] += fading[g]*u  # sum the modulated signal
            end
        else
            for g in 1:G
                u = 2*cword[twoLs+g] - 1  
                signal[g] += u            # sum the modulated signal
            end
        end
    end

end

function resetmatrix!(
    X::Matrix{<:Real},
    Nv::Vector{Vector{Int}},
    value::Real
)    
    
    @inbounds for vj in eachindex(Nv)
        for ci in Nv[vj]
            X[ci,vj] = value
        end
    end
end

# append the CRC to the message
function append_CRC!(
    Cw::Union{Matrix{Bool},Vector{Bool}},
    b::Vector{Bool},
    msg::Vector{Bool},
    g_CRC::Vector{Bool},
    A::Int,
    K::Int
)

    @inbounds begin
        for i in 1:A
            Cw[i] = msg[i]
        end
        for i in A+1:K
            Cw[i] = false
        end
        divide_poly_CRC!(b,Cw,g_CRC,A,K)
    end
end

function print_test(
    text::String,
    array::Vector{Bool}

)    
    println()
    print("$text (L = $(length(array))):")
    for i in eachindex(array)
        if i%80 == 1
            println()
        end
        print(Int(array[i]))
    end
    println()
end 