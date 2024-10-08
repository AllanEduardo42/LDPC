################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Lookup table functions

function lookupTable()

    phi = zeros(SIZE)
    for i in 1:SIZE
        m = exp(i/SIZE_per_RANGE)
        phi[i] = SIZE_per_RANGE*log((m+1)/(m-1))
    end

    return phi
end

function get_index(arg::AbstractFloat)
    
    z = trunc(Int,arg)
    if z >= SIZE
        i = SIZE
    else
        i = z + 1
    end
    
    return i
end