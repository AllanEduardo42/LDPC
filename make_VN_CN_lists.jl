################################################################################
# Allan Eduardo Feitosa
# 26 mai 2025
# Functions to create lists of VNs and CNs

"""For each column j of matrix H, find the indices i where H[i,j] = 1, and 
return a vector "vn2cn" where vn2cn[j] is a vector containing the
indices i where H[i,j] = 1"""
function make_vn2cn_list(H::Matrix{Bool})

    N = size(H,2)
    vn2cn = Vector{Vector{Int}}()
    @inbounds for n in 1:N
        push!(vn2cn,findall(isone, H[:,n]))
    end

    return vn2cn

end

"""For each row i of matrix H, find the indices j where H[i,j] = 1, and 
return a vector "cn2vn" where cn2vn[i] is a vector containing the
indices j where H[i,j] = 1"""
function make_cn2vn_list(H::Matrix{Bool})

    M = size(H,1)
    cn2vn = Vector{Vector{Int}}()
    @inbounds for m in 1:M
        push!(cn2vn,findall(isone, H[m,:]))
    end

    return cn2vn

end