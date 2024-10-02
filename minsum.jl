################################################################################
# Allan Eduardo Feitosa
# 22 set 2024
# Horizontal update of the LLR based MIN SUM Algorithm

function abs_sign!(Lq::AbstractFloat,s::Bool)
    sig = signbit(Lq)
    @fastmath ab = abs(Lq)
    return ab, sig, s ⊻ sig
end

function _minsum!(                           
    Lq::Matrix{<:AbstractFloat},
    signs::Vector{Bool},
    m::Integer,
    cn2vn::Vector{Vector{T}} where {T<:Integer},   
    )

    s = false
    minL = INF
    minL2 = INF
    max_idx = 0
    for n in cn2vn[m]
        @inbounds β, signs[n], s = abs_sign!(Lq[n,m],s)
        if @fastmath β < minL
            max_idx = n
            minL, minL2 = β, minL
        elseif @fastmath β < minL2
            minL2 = β
        end
    end

    return minL, minL2, s, max_idx

end

function 
    __minsum!(
        n::Integer,
        sig::Bool,
        minL::AbstractFloat,
        minL2::AbstractFloat,
        s::Bool,
        max_idx::Integer
    )

    if n == max_idx #(pick the second least Lq)
        return @fastmath (1 - 2*(sig ⊻ s))*minL2
    else
        return @fastmath (1 - 2*(sig ⊻ s))*minL
    end

end

