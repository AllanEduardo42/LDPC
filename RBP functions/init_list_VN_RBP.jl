function
    init_list_VN_RBP!(
        Lq::Matrix{Float64},
        Nc::Vector{Vector{Int}},
        Nv::Vector{Vector{Int}},
        aux::Vector{Float64},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{Float64},Nothing},
        newLr::Matrix{Float64},
        Factors::Vector{Float64},
        inlist::Vector{Bool},
        alpha::Vector{Float64},
        coords::Vector{Int},
        listsize::Int
    )
    
    @inbounds for ci in eachindex(Nc)
        Nci = Nc[ci]
        A, B, C, D = calc_ABCD!(aux,signs,phi,Lq,ci,Nci)
        for vj in Nci
            li = LinearIndices(newLr)[ci,vj]
            newlr = calc_Lr(A,B,C,D,vj,aux,signs,phi)
            newLr[li] = newlr
        end
    end
    @inbounds for vj in eachindex(Nv)
        alp = _calc_all_residues_VN(newLr,vj,Nv[vj])
        alp *= Factors[vj]
        add_list_VN!(alp,alpha,listsize,inlist,coords,vj)
    end
end