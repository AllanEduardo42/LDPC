################################################################################
# Allan Eduardo Feitosa
# 12 May 2025
# Calculate all alpha for the VN-RBP

include("calc_residue.jl")

function
    init_LD_RBP!(
        Lq::Matrix{Float64},
        Nc::Vector{Vector{Int}},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{Float64},Nothing},
        newLr::Matrix{Float64},
        alpha::Vector{Float64},
        Nv::Vector{Vector{Int}}
    )

    @inbounds for ci in eachindex(Nc)
        Nci = Nc[ci]
        for vj in Nci
            newLr[ci,vj] = calc_Lr(Nci,ci,vj,Lq)
        end
    end
    @inbounds for vj in eachindex(Nv)
        alp = 0.0
        for ci in Nv[vj]
            alp += newLr[ci,vj]
        end
        alpha[vj] = abs(alp)
    end
end