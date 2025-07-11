################################################################################
# Allan Eduardo Feitosa
# 26 Feb 2025
# Calculate the residue for the RBP algorithm

function 
    calc_residue!(
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},
        newLr::Matrix{Float64},
        li::Int,
        ci::Int,
        vj::Int,
        Nci::Vector{Int},
        factor::Float64,
        alp::Float64
    )

    oldlr = Lr[li]
    newlr = calc_Lr(Nci,ci,vj,Lq)
    newLr[li] = newlr
    residue = abs(newlr - oldlr)*factor
    if residue > alp
        alp = residue
    end

    return alp, residue
end