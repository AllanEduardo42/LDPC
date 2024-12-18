################################################################################
# Allan Eduardo Feitosa
# 16 dez 2024
# Update check to node messages of the RBP protocol

# 1) if the residues are calculate by TANH or FAST
function 
    RBP_update_Lr!(
        lmax::Integer,
        ::Integer,
        ::Integer,
        ::Vector{Vector{T}} where {T<:Integer},
        ::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Ms::Matrix{<:AbstractFloat},
        ::Vector{<:AbstractFloat},
        ::Nothing,
        ::Nothing
    )

    @inbounds Lr[lmax] = Ms[lmax]

end

function 
    RBP_update_Lr!(
        lmax::Integer,
        ::Integer,
        ::Integer,
        ::Vector{Vector{T}} where {T<:Integer},
        ::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Ms::Matrix{<:AbstractFloat},
        ::Vector{<:AbstractFloat},
        ::Vector{Bool},
        ::Vector{<:AbstractFloat}
    )

    @inbounds Lr[lmax] = Ms[lmax]

end

# 2) otherwise
function 
    RBP_update_Lr!(
        lmax::Integer,
        cnmax::Integer,
        vnmax::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        ::Matrix{<:AbstractFloat},
        ::Union{Vector{<:AbstractFloat},Nothing},
        ::Vector{Bool},
        ::Nothing
    )

    pLr = 1.0
    @fastmath @inbounds for n in cn2vn[cnmax]
        if n != vnmax
            pLr *= tanh(0.5*Lq[n,cnmax])
        end
    end    
    if @fastmath abs(pLr) < 1 
        @fastmath @inbounds Lr[lmax] = 2*atanh(pLr)
    elseif pLr > 0
        @fastmath @inbounds Lr[lmax] = INFFLOAT
    else
        @fastmath @inbounds Lr[lmax] = NINFFLOAT
    end

end