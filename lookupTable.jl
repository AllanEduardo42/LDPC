################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Lookup table functions

function lookupTable()

    phi = zeros(TABLESIZE)
    for i in 1:TABLESIZE
        m = exp(i/SIZE_PER_RANGE)
        phi[i] = SIZE_PER_RANGE*log((m+1)/(m-1))
    end

    return phi
end

function ϕ(Lq::AbstractFloat, phi::Vector{<:AbstractFloat})    
    @inbounds phi[get_index(Lq)]
end

function get_index(arg::AbstractFloat)
    
    if arg > TABLESIZE
        return TABLESIZE
    else
        return trunc(Int,arg) + 1
    end
end