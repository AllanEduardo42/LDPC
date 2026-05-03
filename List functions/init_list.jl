################################################################################
# Allan Eduardo Feitosa
# 3 Mar 2025

include("add_to_list.jl")

function
    init_list!(
        V2C::Matrix{Float64},
        Nc::Vector{Vector{Int}},
        phi::Union{Vector{Float64},Nothing},
        msum_factor::Union{Float64,Nothing},
        newC2V::Matrix{Float64},   
        inlist::Matrix{Bool},
        Residuals::Matrix{Float64},
        list::Vector{Float64},
        coords::Matrix{Int},
        listsize::Int;
        listsize2=listsize
    )

    @fastmath @inbounds for ci in eachindex(Nc)
        Nci = Nc[ci]
        for vj in Nci
            li = LinearIndices(newC2V)[ci,vj]
            newlr = calc_C2V(Nci,ci,vj,V2C,msum_factor)
            newC2V[li] = newlr
            residual = abs(newlr)
            Residuals[li] = residual
            add_to_list!(inlist,list,coords,residual,li,ci,vj,listsize)
        end
    end
end