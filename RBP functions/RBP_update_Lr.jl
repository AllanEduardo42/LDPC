################################################################################
# Allan Eduardo Feitosa
# 23 Feb 2024
# Update check to node messages of the RBP protocol

# 1) if the residues are calculate by the FAST or TABL method
function 
    RBP_update_Lr!(
        lmax::Int,
        Lr::Matrix{Float64},
        newLr::Matrix{Float64},
        ::Int,
        ::Int,
        ::Vector{Int},
        ::Matrix{Float64},        
        ::Vector{Float64},
        ::Union{Vector{Bool},Nothing},
        ::Union{Vector{Float64},Nothing}
    )
    
    # update check to node message Lr[cimax,vjmax]
    @inbounds Lr[lmax] = newLr[lmax]

end

# 2) if the residues are calculate by MSUM method
function 
    RBP_update_Lr!(
        lmax::Int,
        Lr::Matrix{Float64},
        ::Matrix{Float64},
        cimax::Int,
        vjmax::Int,
        Ncimax::Vector{Int},
        Lq::Matrix{Float64},     
        ::Nothing,
        ::Vector{Bool},
        ::Nothing
    )

    # update check to node message Lr[cimax,vjmax]
    @fastmath @inbounds begin
        pLr = 1.0
        for vj in Ncimax
            if vj != vjmax
                lq = Lq[cimax,vj]
                if lq == 0.0
                    return 0.0
                else
                    pLr *= tanh(0.5*lq)
                end
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