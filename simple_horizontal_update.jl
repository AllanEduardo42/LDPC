################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Simple horizontal update of the Sum-Product Algorithm

function simple_horizontal_update(R, M, δq, indices_n)

    for m=1:M
        for n in indices_n[m]
            δr = 1
            for nn in indices_n[m]
                if nn != n
                    @inbounds δr *= δq[nn,m]
                end
            end
            @inbounds R[m,n,1] = 0.5*(1+δr)
            @inbounds R[m,n,2] = 0.5*(1-δr)
        end
    end

    return R

end