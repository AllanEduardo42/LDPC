function
    init_list_VN_RBP!(
        Lq::Matrix{Float64},
        Nc::Vector{Vector{Int}},
        Nv::Vector{Vector{Int}},
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
        for vj in Nci
            li = LinearIndices(newLr)[ci,vj]
            newlr = calc_Lr(Nci,ci,vj,Lq)
            newLr[li] = newlr
        end
    end
    @inbounds for vj in eachindex(Nv)
        alp = 0.0
        for ci in Nv[vj]
            residue = abs(newLr[ci,vj])
            if residue > alp
                alp = residue
            end
        end
        alp *= Factors[vj]
        add_list_VN!(alp,alpha,listsize,inlist,coords,vj)
    end
end