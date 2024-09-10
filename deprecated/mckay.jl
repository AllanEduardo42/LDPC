################################################################################
# Allan Eduardo Feitosa
# 29 ago 2024
# Generate H according to McKay

function mckay(M, N, colord; mincicle=4)

    H = zeros(Int64,M,N)
    density = colord/M
    for n=1:N
        rnd = zeros(M)
        h = zeros(Int64,M)
        while sum(h) != colord
            rand!(rnd)
            h = rnd .< density
        end
        H[:,n] = h
    end   
    
    return H
end
