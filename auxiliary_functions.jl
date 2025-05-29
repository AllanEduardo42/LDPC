################################################################################
# Allan Eduardo Feitosa
# 26 mai 2025
# Auxiliary functions used in sincore.jl

function 
    init_Lq!(
        Lq::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        vn2cn::Vector{Vector{T}} where {T<:Integer}
    )
    
    @inbounds for n in eachindex(vn2cn)
        for m in vn2cn[n]
            Lq[m,n] = Lf[n]
        end
    end
end

function
    init_Lq!(
        Lq::Array{<:AbstractFloat,3},
        Lf::Matrix{<:AbstractFloat},
        vn2cn::Vector{Vector{T}} where {T<:Integer}
    )
   
    @inbounds for n in eachindex(vn2cn)
        for m in vn2cn[n]
            Lq[m,n,1] = Lf[n,1]
            Lq[m,n,2] = Lf[n,2]
        end
    end

end

function
    received_signal!(
        signal::AbstractArray{<:AbstractFloat},
        noise::Vector{<:AbstractFloat},
        G::Integer,
        σ::AbstractFloat,
        rgn_noise::AbstractRNG,
        ::Nothing
    )

    @fastmath begin
        randn!(rgn_noise,noise)
        for i in 1:G
            signal[i] += noise[i]*σ 
        end
    end

end

function
    received_signal!(
        signal::AbstractArray{<:AbstractFloat},
        ::Vector{<:AbstractFloat},
        G::Integer,
        σ::AbstractFloat,
        ::AbstractRNG,
        noisetest::Vector{<:AbstractFloat}
    )

    @fastmath begin
        for i in 1:G
            signal[i] += noisetest[i]*σ
        end
    end

end

function
    resetmatrix!(
        X::Matrix{<:Real},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        value::Real
    )    
    
    @inbounds for n in eachindex(vn2cn)
        for m in vn2cn[n]
            X[m,n] = value
        end
    end
end

function
    resetmatrix!(
        X::Array{<:Real,3},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        value::Real
    )    
    
    @inbounds for n in eachindex(vn2cn)
        for m in vn2cn[n]
            X[m,n,1] = value
            X[m,n,2] = value
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

    msg .= msgtest

end

# append the CRC to the message
function 
    append_CRC!(
        Cw::Union{Matrix{Bool},Vector{Bool}},
        b::Vector{Bool},
        msg::Vector{Bool},
        g_CRC::Vector{Bool},
        A::Integer,
        K::Integer
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