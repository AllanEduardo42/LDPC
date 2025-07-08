################################################################################
# Allan Eduardo Feitosa
# 26 mai 2025
# Auxiliary functions used in sincore.jl

function 
    init_Lq!(
        Lq::Matrix{Float64},
        Lf::Vector{Float64},
        Nv::Vector{Vector{Int}},
        signs::Union{Vector{Bool},Nothing},
        raw::Bool
    )
    
    @inbounds for vj in eachindex(Nv)
        aux = Lf[vj]
        for ci in Nv[vj]
            if raw
                Lq[ci,vj] = aux
            else
                Lq[ci,vj] = tanhLq(aux,signs)
            end
        end
    end
end

function
    init_Lq!(
        Lq::Array{Float64,3},
        Lf::Matrix{Float64},
        Nv::Vector{Vector{Int}}
    )
   
    @inbounds for vj in eachindex(Nv)
        for ci in Nv[vj]
            Lq[ci,vj,1] = Lf[vj,1]
            Lq[ci,vj,2] = Lf[vj,2]
        end
    end

end

function
    received_signal!(
        signal::Vector{Float64},
        cword::Vector{Bool},
        G::Int,
        twoLs::Int,
        σ::Float64,
        rgn::AbstractRNG,
        ::Nothing
    )

    @inbounds @fastmath begin
        randn!(rgn,signal)    # put the noise in the vector 'signal'
        lmul!(σ,signal)             # multiply by the standard deviation
        for g in 1:G
            aux = 2*cword[twoLs+g] - 1  
            signal[g] += aux            # sum the modulated signal
        end
    end

end

function
    received_signal!(
        signal::Vector{Float64},
        cword::Vector{Bool},
        G::Int,
        twoLs::Int,
        σ::Float64,
        ::AbstractRNG,
        noisetest::Vector{Float64}
    )

    @fastmath @fastmath begin
        copy!(signal,noisetest)
        lmul!(σ,signal)
        for g in 1:G
            aux = 2*cword[twoLs+g] - 1
            signal[g] += aux
        end
    end
end

function
    resetmatrix!(
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

function
    resetmatrix!(
        X::Array{<:Real,3},
        Nv::Vector{Vector{Int}},
        value::Real
    )    
    
    @inbounds for vj in eachindex(Nv)
        for ci in Nv[vj]
            X[ci,vj,1] = value
            X[ci,vj,2] = value
        end
    end
end

function generate_message!(
    msg::Vector{Bool},
    rgn_msg::AbstractRNG,
    ::Nothing
)

    rand!(rgn_msg,msg,Bool)

end

function generate_message!(
    msg::Vector{Bool},
    ::AbstractRNG,
    msgtest::Vector{Bool}
)

    copy!(msg, msgtest)

end

# append the CRC to the message
function 
    append_CRC!(
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