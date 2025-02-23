################################################################################
# Allan Eduardo Feitosa
# 23 Feb 2024
# Update check to node messages of the RBP protocol

# 1) if the residues are calculate by TANH and FAST 
function 
    RBP_update_Lr!(
        lmax::Integer,
        Lr::Matrix{<:AbstractFloat},
        Ms::Matrix{<:AbstractFloat},
        ::Integer,
        ::Integer,
        ::Vector{Vector{T}} where {T<:Integer},
        ::Matrix{<:AbstractFloat},        
        ::Union{Vector{<:AbstractFloat},Nothing},
        ::Nothing,
        ::Nothing,
    )
    
    # update check to node message Lr[cnmax,vnmax]
    @inbounds Lr[lmax] = Ms[lmax]

end

# 2) if the residues are calculate by ALTN and TABL
function 
    RBP_update_Lr!(
        lmax::Integer,
        Lr::Matrix{<:AbstractFloat},
        Ms::Matrix{<:AbstractFloat},
        ::Integer,
        ::Integer,
        ::Vector{Vector{T}} where {T<:Integer},
        ::Matrix{<:AbstractFloat},        
        ::Vector{<:AbstractFloat},
        ::Vector{Bool},
        ::Union{Vector{<:AbstractFloat},Nothing},
    )    

    # update check to node message Lr[cnmax,vnmax]
    @inbounds Lr[lmax] = Ms[lmax]

end

# 3) if the residues are calculate by MSUM
function 
    RBP_update_Lr!(
        lmax::Integer,
        Lr::Matrix{<:AbstractFloat},
        ::Matrix{<:AbstractFloat},
        cnmax::Integer,
        vnmax::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        Lq::Matrix{<:AbstractFloat},     
        ::Nothing,
        ::Vector{Bool},
        ::Nothing
    )

    # update check to node message Lr[cnmax,vnmax]
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