################################################################################
# Allan Eduardo Feitosa
# 22 Feb 2024
# Function to find the edge of the maximum residue

function
    findmaxedge(
        Residues::Matrix{Float64},
        alpha::Vector{Float64},
        Nc::Vector{Vector{Int}},
    )

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
                residue = Residues[cimax,vj]
                if residue > maxresidue
                    maxresidue2, maxresidue = maxresidue, residue
                    vjmax = vj
                elseif residue > maxresidue2
                    maxresidue2 = residue
                end
            end
            alpha[cimax] = maxresidue2
            return cimax, vjmax
        end
    end
end

function
    findmaxedge_SVNF(
        Residues::Matrix{Float64},
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
                residue = Residues[ca,vk]
                if residue > maxresidue
                    maxresidue = residue
                    cimax = ca
                    vjmax = vk
                end
            end
        end
    end
    
    return cimax, vjmax

end

function
    findmaxedge_VC(
        Residues::Matrix{Float64},
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
                residue = Residues[ci,vjmax]
                if residue > maxresidue
                    maxresidue2, maxresidue = maxresidue, residue
                    cimax = ci
                elseif residue > maxresidue2
                    maxresidue2 = residue
                end
            end
            alpha[vjmax] = maxresidue2
            return cimax, vjmax
        end
    end
end

