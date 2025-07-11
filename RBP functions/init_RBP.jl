################################################################################
# Allan Eduardo Feitosa
# 3 Mar 2025
# Calculate all Residues

include("calc_residue.jl")

function
    init_RBP!(
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},
        Nc::Vector{Vector{Int}},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{Float64},Nothing},
        newLr::Matrix{Float64}, 
        alpha::Vector{Float64},
        Residues::Matrix{Float64}
    )

    @fastmath @inbounds for ci in eachindex(Nc)
        Nci = Nc[ci]
        alp = 0.0
        for vj in Nci
            li = LinearIndices(Lr)[ci,vj]
            alp, residue = calc_residue!(Lq,Lr,newLr,li,ci,vj,Nci,1.0,alp)
            Residues[li] = residue
        end
        alpha[ci] = alp
    end
end