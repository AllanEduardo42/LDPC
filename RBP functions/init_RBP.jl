################################################################################
# Allan Eduardo Feitosa
# 3 Mar 2025
# Calculate all Residues

include("calc_residue.jl")

# FAST, TABL and MSUM
function
    init_RBP!(
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},
        Nc::Vector{Vector{Int}},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{Float64},Nothing},
        newLr::Matrix{Float64},
        Factors::Matrix{Float64},        
        alpha::Vector{Float64},
        Residues::Matrix{Float64},
        relative::Bool,
        raw::Bool
    )

    @fastmath @inbounds if raw
        for ci in eachindex(Nc)
            Nci = Nc[ci]
            alp = 0.0
            for vj in Nci
                li = LinearIndices(newLr)[ci,vj]
                newlr = calc_Lr(Nci,ci,vj,Lq)
                newLr[li] = newlr
                residue = abs(newlr - Lr[li])
                Residues[li] = residue
                if residue > alp
                    alp = residue
                end
            end
            alpha[ci] = alp
        end
    else    
        for ci in eachindex(Nc)
            Nci = Nc[ci]
            A, B, C, D = calc_ABCD!(Lq,ci,Nci,signs,phi)
            alp = 0.0
            for vj in Nci
                li = LinearIndices(newLr)[ci,vj]
                newlr = calc_Lr(A,B,C,D,vj,Lq[li],signs,phi)
                newLr[li] = newlr
                residue = calc_residue(newLr[li],Lr[li],Factors[li],relative,Lq[li])
                Residues[li] = residue
                if residue > alp
                    alp = residue
                end
            end
            alpha[ci] = alp
        end
    end
end