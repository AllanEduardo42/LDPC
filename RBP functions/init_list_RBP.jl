################################################################################
# Allan Eduardo Feitosa
# 3 Mar 2025
# Calculate all residues

include("calc_residue.jl")
include("add_residue.jl")

# FAST, TABL and MSUM
function
    init_list_RBP!(
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},
        Nc::Vector{Vector{Int}},
        aux::Vector{Float64},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{Float64},Nothing},
        newLr::Matrix{Float64},
        Factors::Matrix{Float64},        
        inlist::Union{Matrix{Int},Matrix{Bool}},
        residues::Vector{Float64},
        coords::Matrix{Int},
        listsizes::Vector{Int}
    )
    
    @fastmath @inbounds for ci in eachindex(Nc)
        Nci = Nc[ci]
        A, B, C, D = calc_ABCD!(aux,signs,phi,Lq,ci,Nci)
        for vj in Nci
            li = LinearIndices(newLr)[ci,vj]
            newlr = calc_Lr(A,B,C,D,vj,aux,signs,phi)
            newLr[li] = newlr
            residue = calc_residue(newLr[li],Lr[li],Factors[li],false,Lq[li])
            add_residue!(inlist,residues,coords,residue,li,ci,vj,listsizes[1])
        end
    end
end

# TANH
function 
    init_list_RBP!(
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},
        Nc::Vector{Vector{Int}},
        ::Nothing,
        ::Nothing,
        ::Nothing,
        newLr::Matrix{Float64},
        Factors::Matrix{Float64},        
        inlist::Union{Matrix{Int},Matrix{Bool}},
        residues::Vector{Float64},
        coords::Matrix{Int},
        listsizes::Vector{Int}
    )

    @fastmath @inbounds for ci in eachindex(Nc)
        Nci = Nc[ci]
        for vj in Nci
            li = LinearIndices(newLr)[ci,vj]
            newlr = calc_Lr(Nci,ci,vj,Lq)
            newLr[li] = newlr
            residue = calc_residue_raw(newLr[li],Lr[li],Factors[li],false,Lq[li])
            add_residue!(inlist,residues,coords,residue,li,ci,vj,listsizes[1])
        end
    end
end