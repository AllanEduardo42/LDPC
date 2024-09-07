################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Lookup table functions

function lookupTable()
    phi = zeros(SIZE)
    for i=1:SIZE
        m = exp(i/SIZE_per_RANGE)
        phi[i] = SIZE_per_RANGE*log((m+1)/(m-1))
    end

    return phi
end

function get_phi(arg::Float64, phi::Vector{Float64})::Float64
    z = unsafe_trunc(Int,arg)
    if z >= SIZE
        return 0.0
    else
        return phi[z + 1]
    end    
end