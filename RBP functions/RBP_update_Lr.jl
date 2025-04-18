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
        ::Vector{<:Integer},
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
        ::Vector{<:Integer},
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
        vns::Vector{<:Integer},
        Lq::Matrix{<:AbstractFloat},     
        ::Nothing,
        ::Vector{Bool},
        ::Nothing
    )

    # update check to node message Lr[cnmax,vnmax]
    @fastmath @inbounds begin
        pLr = 1.0
        for n in vns
            if n != vnmax
                pLr *= tanh(0.5*Lq[n,cnmax])
            end
        end    
        if abs(pLr) < 1 
            Lr[lmax] = 2*atanh(pLr)
        elseif pLr > 0
            Lr[lmax] = INFFLOAT
        else
            Lr[lmax] = NINFFLOAT
        end
    end

end