################################################################################
# Allan Eduardo Feitosa
# 22 Feb 2024
# auxiliary functions for RBP-based algorithm

# Calculates all the residues at the algorithm initialization
function init_residuals!(
    V2C::Matrix{Float64},
    Nc::Vector{Vector{Int}},
    newC2V::Matrix{Float64}, 
    alpha::Vector{Float64},
    Residuals::Matrix{Float64},
    msum_factor::Union{Float64,Nothing}
)

    # for ci in eachindex(Nc)
    @fastmath @inbounds for ci in eachindex(Nc)
        Nci = Nc[ci]
        alp = 0.0
        for vj in Nci
            li = LinearIndices(newC2V)[ci,vj]
            newc2v = calc_C2V(Nci,ci,vj,V2C,msum_factor)
            newC2V[li] = newc2v
            residual = abs(newc2v)
            if residual > alp
                alp = residual
            end
            Residuals[li] = residual
        end
        alpha[ci] = alp
    end
end


# Function to find the edge of the maximum residual
function findmaxedge(
    Residuals::Matrix{Float64},
    alpha::Vector{Float64},
    Nc::Vector{Vector{Int}}
)

    # begin
    @fastmath @inbounds begin
        cimax = 0
        maxalp = 0.0
        for ci in eachindex(alpha)
            alp = alpha[ci]
            if alp > maxalp
                maxalp = alp
                cimax = ci
            end
        end
        
        if cimax == 0
            return 0, 0
        else
            vjmax = 0
            maxresidue = 0.0
            maxresidue2 = 0.0
            for vj in Nc[cimax]            
                residual = Residuals[cimax,vj]
                if residual > maxresidue
                    maxresidue2, maxresidue = maxresidue, residual
                    vjmax = vj
                elseif residual > maxresidue2
                    maxresidue2 = residual
                end
            end
            alpha[cimax] = maxresidue2
            return cimax, vjmax
        end
    end
end

# Function to find the edge of the maximum residual for sVNF
function findmaxedge_SVNF(
    Residuals::Matrix{Float64},
    vj::Integer,
    Nvj::Vector{Int},
    Nc::Vector{Vector{Int}}
)

    cimax = 0
    vjmax = 0
    maxresidue = 0.0
    @fastmath @inbounds for ca in Nvj
        Nca = Nc[ca]
        for vk in Nca
            if vk ≠ vj
                residual = Residuals[ca,vk]
                if residual > maxresidue
                    maxresidue = residual
                    cimax = ca
                    vjmax = vk
                end
            end
        end
    end
    
    return cimax, vjmax

end

# Function to find the edge of the maximum residual for UBP
function findmaxedge_UBP(
    Residuals::Matrix{Float64},
    alpha::Vector{Float64},
    Nc::Vector{Vector{Int}},
    UBP::Vector{Bool}
)

    # begin
    @fastmath @inbounds begin
        cimax = 0
        maxalp = 0.0
        for ci in eachindex(alpha)
            if UBP[ci] == true
                alp = alpha[ci]
                if alp > maxalp
                    maxalp = alp
                    cimax = ci
                end
            end
        end
        
        if cimax == 0
            return 0, 0
        else
            vjmax = 0
            maxresidue = 0.0
            maxresidue2 = 0.0
            for vj in Nc[cimax]            
                residual = Residuals[cimax,vj]
                if residual > maxresidue
                    maxresidue2, maxresidue = maxresidue, residual
                    vjmax = vj
                elseif residual > maxresidue2
                    maxresidue2 = residual
                end
            end
            alpha[cimax] = maxresidue2
            return cimax, vjmax
        end
    end
end

# Function to find the edge of the maximum residual for RBP-D1VN
function findmaxedge_D1VN(
    Residuals::Matrix{Float64},
    alpha::Vector{Float64},
    Nc::Vector{Vector{Int}},
    F::Vector{Bool}
)

    # begin
    @fastmath @inbounds begin
        maxresidue = 0.0
        cimax = 0
        vjmax = 0
        for ci in eachindex(Nc)
            Nci = Nc[ci]
            for vj in Nci
                if !F[vj]
                    residual = Residuals[ci,vj]
                    if residual > maxresidue
                        maxresidue = residual
                        cimax = ci
                        vjmax = vj
                    end
                end
            end
        end
    end

    return cimax, vjmax
end

# Function to find the check node with the maximum residual for NW-RBP
function findmaxnode(
    alpha::Vector{Float64}        
)

    cimax = 0
    max_alp = 0.0
    @inbounds @fastmath for e in eachindex(alpha)
        alp = alpha[e]
        if alp > max_alp
            max_alp = alp
            cimax = e
        end
    end

    return cimax

end
