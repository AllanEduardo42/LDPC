function llr_simple_horizontal_update(Lr::Matrix{Float64}, M, Lq::Matrix{Float64},
                                      indices_n::Vector{Vector{Int64}})

    for m=1:M
        for n in indices_n[m]
            pLr = 1.0
            for nn in indices_n[m]
                if nn != n
                    pLr *= tanh(0.5*Lq[nn,m])
                end
            end
            Lr[m,n] = 2*atanh(pLr)
        end
    end

    return Lr

end

function llr_simple_horizontal_update_alt(Lr::Matrix{Float64}, M, Lq::Matrix{Float64},
                                          indices_n::Vector{Vector{Int64}})

    for m=1:M
        for n in indices_n[m]
            sum = 0.0
            s = 1
            for nn in indices_n[m]
                if nn != n
                    sum += ϕ(abs(Lq[nn,m]))
                    s *= sign(Lq[nn,m])
                end
            end
            Lr[m,n] = s*ϕ(abs(sum))
        end
    end

    return Lr

end

function ϕ(x)
    m = exp(x)
    return log(m + 1) - log(m - 1)
end