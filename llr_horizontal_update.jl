################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Horizontal update of the LLR based Sum-Product Algorithm

function llr_horizontal_update(Lr::Matrix{Float64}, M::Int64,
                                      Lq::Matrix{Float64},
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

# function llr_horizontal_update_alt(Lr::Matrix{Float64}, M::Int64,
#                                           Lq::Matrix{Float64},
#                                           indices_n::Vector{Vector{Int64}},
#                                           Lrn::Vector{Float64},
#                                           sn::Vector{Int64})
#     for m=1:M
#         sum = 0.0
#         s = 1 
#         for n in indices_n[m]
#             Lrn[n] = ϕ(abs(Lq[n,m]))
#             sum += Lrn[n]
#             sn[n] = sign(Lq[n,m])
#             s *= sn[n]
#         end
#         for n in indices_n[m]
#             Lr[m,n] = (s*sn[n])*(ϕ(abs(sum - Lrn[n])))
#         end    
#     end

#     return Lr

# end

function llr_horizontal_update_alt(Lr::Matrix{Float64}, M::Int64,
                                         Lq::Matrix{Float64},
                                         indices_n::Vector{Vector{Int64}},
                                         Lrn::Vector{Float64},
                                         sn::Vector{Int64})

    for m=1:M
        sum = 0.0
        s = 1 
        for n in indices_n[m]
            if Lq[n,m] >= 0
                Lrn[n] = ϕ(Lq[n,m])
                sn[n] = 1
            else
                Lrn[n] = ϕ(-Lq[n,m])
                sn[n] = -1
                s *= sn[n]
            end
            sum += Lrn[n]
        end
        for n in indices_n[m]
            if sn[n] != s
                Lr[m,n] = -ϕ(abs(sum - Lrn[n]))
            else
                Lr[m,n] = ϕ(abs(sum - Lrn[n]))
            end
        end    
    end

    return Lr

end

function llr_horizontal_update_table(Lr::Matrix{Float64}, M::Int64,
                                            Lq::Matrix{Float64},
                                            indices_n::Vector{Vector{Int64}},
                                            Lrn::Vector{Float64},
                                            sn::Vector{Int64},
                                            phi::Vector{Float64})
    for m=1:M
        sum = 0.0
        s = 1 
        for n in indices_n[m]
            if Lq[n,m] >= 0
                Lrn[n] = phi[get_index(Lq[n,m])]
                sn[n] = 1
            else
                Lrn[n] = phi[get_index(-Lq[n,m])]
                sn[n] = -1
                s *= sn[n]
            end
            sum += Lrn[n]
        end
        for n in indices_n[m]
            if sn[n] != s
                Lr[m,n] = -phi[get_index(abs(sum - Lrn[n]))]
            else
                Lr[m,n] = phi[get_index(abs(sum - Lrn[n]))]
            end
        end    
    end

    return Lr

end

function ϕ(x)
    m = exp(x)-1
    return log(1 + 2/m)
end

function get_index(arg::Float64)::Int64
    z = unsafe_trunc(Int,arg)
    if z >= SIZE
        i = SIZE
    else
        i = z + 1
    end
    
    return i
end