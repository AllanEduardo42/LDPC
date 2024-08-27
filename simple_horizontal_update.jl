################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Simple horizontal update of the Sum-Product Algorithm

function simple_horizontal_update(R, M, q, indices_n)

    δq = q[:,:,1]-q[:,:,2]

    for m=1:M
        for n in indices_n[m]
            δr = 1
            for nn in indices_n[m]
                if nn != n
                    δr *= δq[nn,m]
                end
            end
            R[m,n,1] = 0.5*(1+δr)
            R[m,n,2] = 0.5*(1-δr)
        end
    end

    return R

end