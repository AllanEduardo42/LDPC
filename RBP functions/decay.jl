################################################################################
# Allan Eduardo Feitosa
# 23 Feb 2024
# function to decay the residuals

function decay!(
    cnmax::Integer,
    vnmax::Integer,
    Factors::Matrix{<:AbstractFloat},
    decayfactor::AbstractFloat,
)

    # decay the RBP factor corresponding to the maximum residue
    @inbounds lmax = LinearIndices(Factors)[cnmax,vnmax] # for optimization
    @fastmath @inbounds Factors[lmax] *= decayfactor

    return lmax
end