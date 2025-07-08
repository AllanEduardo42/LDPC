################################################################################
# Allan Eduardo Feitosa
# 12 May 2025
# Calculate all alpha for the VN-RBP

include("calc_residue.jl")

function
    init_VN_RBP!(
        Lq::Matrix{Float64},
        Nc::Vector{Vector{Int}},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{Float64},Nothing},
        newLr::Matrix{Float64},
        alpha::Vector{Float64},
        Nv::Vector{Vector{Int}},
        raw::Bool
    )

    @inbounds if raw
        for ci in eachindex(Nc)
            Nci = Nc[ci]
            for vj in Nci
                li = LinearIndices(newLr)[ci,vj]
                newlr = calc_Lr(Nci,ci,vj,Lq)
                newLr[li] = newlr
            end
        end
    else    
        for ci in eachindex(Nc)
            Nci = Nc[ci]
            A, B, C, D = calc_ABCD!(Lq,ci,Nci,signs,phi)
            for vj in Nci
                li = LinearIndices(newLr)[ci,vj]
                newlr = calc_Lr(A,B,C,D,vj,Lq[li],signs,phi)
                newLr[li] = newlr
            end
        end
    end
    @inbounds for vj in eachindex(Nv)
        alpha[vj] = _calc_all_residues_VN(newLr,vj,Nv[vj])
    end
end

# RAW
function
    init_VN_RBP!(
        Lq::Matrix{Float64},
        Nc::Vector{Vector{Int}},
        ::Nothing,
        ::Nothing,
        ::Nothing,
        Lr::Matrix{Float64},
        newLr::Matrix{Float64},
        alpha::Vector{Float64},
        Nv::Vector{Vector{Int}}
    )
    
    @inbounds for ci in eachindex(Nc)
        Nci = Nc[ci]
        for vj in Nci
            li = LinearIndices(Lr)[ci,vj]
            newlr = calc_Lr(Nci,ci,vj,Lq)
            newLr[li] = newlr
            # Lr[li] = newlr
        end
    end
    @inbounds for vj in eachindex(Nv)
        alpha[vj] = _calc_all_residues_VN(newLr,vj,Nv[vj])
    end
end

function 
    _calc_all_residues_VN(
        newLr::Matrix{Float64},
        vj::Int,
        Nvj::Vector{Int}
    )

    residue = 0.0   
    @inbounds @fastmath for ci in Nvj
        aux = abs(newLr[ci,vj])
        if  aux > residue
            residue = aux
        end
    end

    return abs(residue)
end