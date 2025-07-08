################################################################################
# Allan Eduardo Feitosa
# 12 May 2025
# Calculate all alpha for the NW-RBP

include("calc_residue.jl")

function
    init_NW_RBP!(
        Lq::Matrix{Float64},
        Nc::Vector{Vector{Int}},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{Float64},Nothing},
        newLr::Matrix{Float64},
        alpha::Vector{Float64},
        raw::Bool
    )

    @fastmath @inbounds if raw
        for ci in eachindex(Nc)
            Nci = Nc[ci]
            maxresidue = 0.0
            for vj in Nci
                li = LinearIndices(newLr)[ci,vj]
                newlr = calc_Lr(Nci,ci,vj,Lq)
                newLr[li] = newlr
                residue = calc_residue_VN_NW_raw(newlr,0.0,1.0)
                if residue > maxresidue
                    maxresidue = residue
                end
            end
            alpha[ci] = maxresidue
        end
    else    
        for ci in eachindex(Nc)
            Nci = Nc[ci]
            A, B, C, D = calc_ABCD!(Lq,ci,Nci,signs,phi)
            maxresidue = 0.0
            for vj in Nci
                li = LinearIndices(newLr)[ci,vj]
                newlr = calc_Lr(A,B,C,D,vj,Lq[li],signs,phi)
                newLr[li] = newlr
                residue = abs(newlr)
                if residue > maxresidue
                    maxresidue = residue
                end
            end
            alpha[ci] = maxresidue
        end
    end
end