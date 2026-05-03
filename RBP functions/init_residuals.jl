################################################################################
# Allan Eduardo Feitosa
# 3 Mar 2025
# Calculate all Residuals

function
    init_residuals!(
        V2C::Matrix{Float64},
        Nc::Vector{Vector{Int}},
        phi::Union{Vector{Float64},Nothing},
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
            residue = abs(newc2v)
            if residue > alp
                alp = residue
            end
            Residuals[li] = residue
        end
        alpha[ci] = alp
    end
end