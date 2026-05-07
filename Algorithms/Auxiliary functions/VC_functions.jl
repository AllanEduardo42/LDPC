################################################################################
# Allan Eduardo Feitosa
# 07 Mai 2026
# auxiliary functions for VC_RBP

function init_VC!(
    V2C::Matrix{Float64},
    Nv::Vector{Vector{Int}},
    alpha::Vector{Float64},
    Residuals::Matrix{Float64}
)

    @inbounds @fastmath for vj in eachindex(Nv)
        alp = 0.0
        for ci in Nv[vj]
            li = LinearIndices(V2C)[ci,vj]
            residual = abs(V2C[li])
            if residual > alp
                alp = residual
            end
            Residuals[li] = residual
        end
        alpha[vj] = alp
    end

end

# Function to find the edge of the maximum residual for VC
function findmaxedge_VC(
    Residuals::Matrix{Float64},
    alpha::Vector{Float64},
    Nv::Vector{Vector{Int}},
)

    @fastmath @inbounds begin
        vjmax = 0
        maxalp = 0.0
        for vj in eachindex(alpha)
            alp = alpha[vj]
            if alp > maxalp
                maxalp = alp
                vjmax = vj
            end
        end
        
        if vjmax == 0
            return 0, 0
        else
            cimax = 0
            maxresidue = 0.0
            maxresidue2 = 0.0
            for ci in Nv[vjmax]            
                residual = Residuals[ci,vjmax]
                if residual > maxresidue
                    maxresidue2, maxresidue = maxresidue, residual
                    cimax = ci
                elseif residual > maxresidue2
                    maxresidue2 = residual
                end
            end
            alpha[vjmax] = maxresidue2
            return cimax, vjmax
        end
    end
end