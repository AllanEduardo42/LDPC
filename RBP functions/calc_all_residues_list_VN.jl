include("add_residue_VN.jl")


function
    calc_all_residues_list_VN!(
        Lq::Matrix{<:AbstractFloat},
        Nc::Vector{Vector{T}} where {T<:Integer},
        aux::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{<:AbstractFloat},Nothing},
        Lr::Matrix{<:AbstractFloat},
        newLr::Matrix{<:AbstractFloat},
        alpha::Vector{<:AbstractFloat},
        Nv::Vector{Vector{T}} where {T<:Integer},
        inlist::Vector{Bool},
        listsizes::Vector{<:Integer},
        coords::Vector{<:Integer}
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
         residue = _calc_all_residues_VN(newLr,vj,Nv[vj])
         add_residue_VN!(inlist,alpha,coords,residue,vj,listsizes[1])
    end
end