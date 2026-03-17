################################################################################
# Allan Eduardo Feitosa
# 12 May 2025
# Calculate all alpha for the NW-RBP

include("calc_residue.jl")

function
    init_NW_RBP!(
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},
        Nc::Vector{Vector{Int}},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{Float64},Nothing},
        newLr::Matrix{Float64},
        alpha::Vector{Float64}
    )

    @fastmath @inbounds for ci in eachindex(Nc)
        Nci = Nc[ci]
        alp = 0.0
        for vj in Nci
            li = LinearIndices(Lq)[ci,vj]
            alp,_ = calc_residue!(Lq,Lr,newLr,li,ci,vj,Nci,1.0,alp)
        end
        alpha[ci] = alp
    end
end