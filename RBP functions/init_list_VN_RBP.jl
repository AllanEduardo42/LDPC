function
    init_list_VN_RBP!(
        Lq::Matrix{Float64},
        Nc::Vector{Vector{Int}},
        aux::Vector{Float64},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{Float64},Nothing},
        Lr::Matrix{Float64},
        newLr::Matrix{Float64},
        alpha::Vector{Float64},
        Nv::Vector{Vector{Int}},
        inlist::Vector{Bool},
        listsizes::Vector{Int},
        coords::Vector{Int}
    )
    
    @inbounds begin
        for ci in eachindex(Nc)
            Nci = Nc[ci]
            A, B, C, D = calc_ABCD!(aux,signs,phi,Lq,ci,Nci)
            for vj in Nci
                li = LinearIndices(Lr)[ci,vj]
                newlr = calc_Lr(A,B,C,D,vj,aux,signs,phi)
                newLr[li] = newlr
                Lr[li] = newlr
            end
        end
        for vj in eachindex(Nv)
            residue = _calc_all_residues_VN(newLr,vj,Nv[vj],false)
            add_residue_VN!(inlist,alpha,coords,residue,vj,listsizes[1])
        end
    end
end