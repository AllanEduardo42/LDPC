################################################################################
# Allan Eduardo Feitosa
# 22 set 2024
# Horizontal update of the LLR based MIN SUM Algorithm

function abs_sign!(Lq::AbstractFloat,s::Bool)
    signs = signbit(Lq)
    return abs(Lq), signs, s ⊻ signs
end

function _minsum!(                           
    Lq::Matrix{<:AbstractFloat},
    signs::Vector{Bool},
    m::Integer,
    cn2vn::Vector{Vector{T}} where {T<:Integer},   
    )

    s = false
    minL = Inf
    minL2 = Inf
    max_idx = 0
    for n in cn2vn[m]
        @inbounds β, signs[n], s = abs_sign!(Lq[n,m],s)
        if β < minL
            max_idx = n
            minL, minL2 = β, minL
        elseif β < minL2
            minL2 = β
        end
    end

    return minL, minL2, s, max_idx

end

function 
    __minsum!(
        n::Integer,
        signs::Bool,
        minL::AbstractFloat,
        minL2::AbstractFloat,
        s::Bool,
        max_idx::Integer
    )

    if n == max_idx #(pick the second least Lq)
        return (1 - 2*(signs ⊻ s))*minL2
    else
        return (1 - 2*(signs ⊻ s))*minL
    end

end

