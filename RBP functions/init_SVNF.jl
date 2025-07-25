################################################################################
# Allan Eduardo Feitosa
# 27 Jun 2025
# Calculate all residues

include("calc_residue.jl")

# FAST, TABL and MSUM
function
    init_SVNF!(
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},
        Nc::Vector{Vector{Int}},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{Float64},Nothing},
        newLr::Matrix{Float64},
        residues::Matrix{Float64}
    )
    
    @fastmath @inbounds for ci in eachindex(Nc)
        Nci = Nc[ci]
        for vj in Nci
            li = LinearIndices(newLr)[ci,vj]
            newlr = calc_Lr(Nci,ci,vj,Lq)
            newLr[li] = newlr
            residues[li] = abs(newlr - Lr[li])
        end
    end
end