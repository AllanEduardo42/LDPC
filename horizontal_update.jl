################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Horizontal update of the Sum-Product Algorithm

function horizontal_update(M, N, q, indices_n)

    r = zeros(M,N,2)

    for m=1:M
        S = length(indices_n[m])-1
        for n in indices_n[m]            
            for s = 0:2^S-1
                dig = digits(s, base = 2, pad = S)
                count = 0
                rr = 1
                if iseven(sum(dig))
                    for nn in indices_n[m] 
                        if nn != n
                            count += 1
                            rr *= q[nn,m,dig[count]+1]
                        end
                    end
                    r[m,n,1] += rr
                else
                    for nn in indices_n[m]
                        if nn != n
                            count += 1
                            rr *= q[nn,m,dig[count]+1]
                        end               
                    end
                    r[m,n,2] += rr
                end
            end
        end
    end

    return r

end