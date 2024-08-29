################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Horizontal update of the LLR based Sum-Product Algorithm

function llr_horizontal_update(M::Int64,
                               Lr::Matrix{Float64},
                               Lq::Matrix{Float64},
                               indices_n::Vector{Vector{Int64}})

@inbounds for m=1:M
        @inbounds for n in indices_n[m]
            pLr = 1.0
            @inbounds for nn in indices_n[m]
                if nn != n
                    @fastmath pLr *= tanh(0.5*Lq[nn,m])
                end
            end
            @fastmath Lr[m,n] = 2*atanh(pLr)
        end
    end

    return Lr

end

function llr_horizontal_update(M::Int64,
                               Lr::Matrix{Float64},                               
                               Lq::Matrix{Float64},
                               indices_n::Vector{Vector{Int64}},
                               sn::Vector{Int64})    

    @inbounds for m=1:M
        minL = Inf
        minL2 = Inf
        idx = 0
        s = 1
        @inbounds for n in indices_n[m]
            Lrn = Lq[n,m]
            if Lrn > 0
                sn[n] = 1
                if Lrn < minL
                    idx = n
                    minL2 = minL
                    minL = Lrn
                elseif Lrn < minL2
                    minL2 = Lrn
                end
            else
                sn[n] = -1
                s *= sn[n]
                Lrn *= -1
                if Lrn < minL
                    idx = n
                    minL2 = minL
                    minL = Lrn
                elseif Lrn < minL2
                    minL2 = Lrn
                end
            end
        end
        @inbounds for n in indices_n[m]
            if n == idx
                if sn[n] != s
                    Lr[m,n] = -minL2
                else
                    Lr[m,n] = minL2
                end
            else
                if sn[n] != s
                    Lr[m,n] = -minL
                else
                    Lr[m,n] = minL
                end
            end
        end  
    end

    return Lr

end

function llr_horizontal_update(M::Int64,
                               Lr::Matrix{Float64},
                               Lq::Matrix{Float64},
                               indices_n::Vector{Vector{Int64}},
                               sn::Vector{Int64},
                               Lrn::Vector{Float64})

    @inbounds for m=1:M
        sum = 0.0
        s = 1 
        @inbounds for n in indices_n[m]
            ####### what it does: #######
            # Lrn[n] = ϕ(abs(Lq[n,m]))
            # sum += Lrn[n]
            # sn[n] = sign(Lq[n,m])
            # s *= sn[n]
            #############################
            if Lq[n,m] >= 0
                Lrn[n] = ϕ(Lq[n,m])
                sn[n] = 1
            else
                Lrn[n] = ϕ(-Lq[n,m])
                sn[n] = -1
                s *= sn[n]
            end
            @fastmath sum += Lrn[n]
        end
        @inbounds for n in indices_n[m]
            ############### what it does ###############
            # Lr[m,n] = (s*sn[n])*(ϕ(abs(sum - Lrn[n])))
            ############################################
            if sn[n] != s
                @fastmath Lr[m,n] = -ϕ(abs(sum - Lrn[n]))
            else
                @fastmath Lr[m,n] = ϕ(abs(sum - Lrn[n]))
            end
        end    
    end

    return Lr

end

function llr_horizontal_update(M::Int64,
                               Lr::Matrix{Float64},
                               Lq::Matrix{Float64},
                               indices_n::Vector{Vector{Int64}},
                               sn::Vector{Int64},
                               Lrn::Vector{Float64},
                               phi::Vector{Float64})
    @inbounds for m=1:M
        sum = 0.0
        s = 1 
        @inbounds for n in indices_n[m]
            if Lq[n,m] >= 0
                Lrn[n] = phi[get_index(Lq[n,m])]
                sn[n] = 1
            else
                Lrn[n] = phi[get_index(-Lq[n,m])]
                sn[n] = -1
                s *= sn[n]
            end
            @fastmath sum += Lrn[n]
        end
        @inbounds for n in indices_n[m]
            if sn[n] != s
                @fastmath Lr[m,n] = -phi[get_index(abs(sum - Lrn[n]))]
            else
                @fastmath Lr[m,n] = phi[get_index(abs(sum - Lrn[n]))]
            end
        end    
    end

    return Lr

end


function ϕ(x::Float64)::Float64
    @fastmath m = exp(x)-1
    @fastmath log(1 + 2/m)
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

