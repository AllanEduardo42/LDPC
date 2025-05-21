################################################################################
# Allan Eduardo Feitosa
# 12 May 2025
# Calculate all alpha for the VN-RBP

include("calc_residue.jl")

function
    calc_all_residues_VN!(
        Lq::Matrix{<:AbstractFloat},
        Nc::Vector{Vector{T}} where {T<:Integer},
        aux::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{<:AbstractFloat},Nothing},
        Lr::Matrix{<:AbstractFloat},
        newLr::Matrix{<:AbstractFloat},
        alpha::Vector{<:AbstractFloat},
        Nv::Vector{Vector{T}} where {T<:Integer},
    )
    
    @inbounds for ci in eachindex(Nc)
        Nci = Nc[ci]
        A, B, C, D = calc_ABCD!(aux,signs,phi,Lq,ci,Nci)
        for vj in Nci
            li = LinearIndices(Lr)[ci,vj]
            newlr = calc_Lr(A,B,C,D,vj,aux,signs,phi)
            newLr[li] = newlr
            Lr[li] = newlr
        end
    end
    @inbounds for vj in eachindex(Nv)
        alpha[vj] = _calc_all_residues_VN(newLr,vj,Nv[vj])
    end
end

# RAW
function
    calc_all_residues_VN!(
        Lq::Matrix{<:AbstractFloat},
        Nc::Vector{Vector{T}} where {T<:Integer},
        ::Nothing,
        ::Nothing,
        ::Nothing,
        Lr::Matrix{<:AbstractFloat},
        newLr::Matrix{<:AbstractFloat},
        alpha::Vector{<:AbstractFloat},
        Nv::Vector{Vector{T}} where {T<:Integer},
    )
    
    @inbounds for ci in eachindex(Nc)
        Nci = Nc[ci]
        for vj in Nci
            li = LinearIndices(Lr)[ci,vj]
            newlr = calc_Lr(Nci,ci,vj,Lq)
            newLr[li] = newlr
            Lr[li] = newlr
        end
    end
    @inbounds for vj in eachindex(Nv)
        alpha[vj] = _calc_all_residues_VN(newLr,vj,Nv[vj])
    end
end

function 
    _calc_all_residues_VN(
        newLr::Matrix{<:AbstractFloat},
        vj::Integer,
        Nvj::Vector{<:Integer}
    )

    residue = 0.0
    for ci in Nvj          
        residue += newLr[ci,vj]
    end

    return abs(residue)
end