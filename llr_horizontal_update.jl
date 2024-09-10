################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Horizontal update of the LLR based Sum-Product Algorithm

using SparseArrays

# TANH
function 
    llr_horizontal_update!(
        Lr::Matrix{Float64},
        Lq::Matrix{Float64},
        indices_row::Vector{Vector{Int64}},
        x::Nothing,
        y::Nothing,
        z::Nothing
    )

    m = 0
    for indices in indices_row
        m += 1
        for n in indices
            pLr = 1.0
            for nn in indices
                if nn != n
                    @inbounds pLr *= tanh(0.5*Lq[nn,m])
                end
            end
            @inbounds Lr[m,n] = 2*atanh(pLr)
        end
    end
end

# APPROX
function
    llr_horizontal_update!(
        Lr::Matrix{Float64},                           
        Lq::Matrix{Float64},
        indices_row::Vector{Vector{Int64}},
        sn::Vector{Bool},
        x::Nothing,
        y::Nothing
    )    

    m = 0
    for indices in indices_row
        m += 1
        s = true
        minL = 0.0
        minL2 = 0.0
        idx = 0
        for n in indices
            @inbounds Lrn = Lq[n,m]
            if Lrn > 0
                @inbounds sn[n] = false
            else
                @inbounds sn[n] = true
                @inbounds s ⊻= sn[n]
                @fastmath Lrn *= -1
            end
            if minL == 0.0
                idx = n
                minL = Lrn
            elseif minL2 == 0
                if Lrn < minL
                    idx = n
                    minL2 = minL
                    minL = Lrn
                else
                    minL2 = Lrn
                end
            else
                if Lrn < minL
                    idx = n
                    minL2 = minL
                    minL = Lrn
                elseif Lrn < minL2
                    minL2 = Lrn
                end
            end
        end
        for n in indices
            if n == idx
                if @inbounds sn[n] == s
                    @inbounds Lr[m,n] = -minL2
                else
                    @inbounds Lr[m,n] = minL2
                end
            else
                if @inbounds sn[n] == s
                    @inbounds Lr[m,n] = -minL
                else
                    @inbounds Lr[m,n] = minL
                end
            end
        end  
    end
end

# ALT
function 
    llr_horizontal_update!(
        Lr::Matrix{Float64},                            
        Lq::Matrix{Float64},
        indices_row::Vector{Vector{Int64}},
        sn::Vector{Bool},
        Lrn::Vector{Float64},
        x::Nothing
    )

    m = 0
    for indices in indices_row
        m += 1
        sum = 0.0
        s = 1 
        for n in indices
            ####### what it does: #######
            # Lrn[n] = ϕ(abs(Lq[n,m]))
            # sum += Lrn[n]
            # sn[n] = sign(Lq[n,m])
            # s *= sn[n]
            #############################
            if @inbounds Lq[n,m] >= 0
                @inbounds Lrn[n] = ϕ(Lq[n,m])
                @inbounds sn[n] = false
            else
                @inbounds Lrn[n] = ϕ(-Lq[n,m])
                @inbounds sn[n] = true
                @inbounds s ⊻= sn[n]
            end
            sum += Lrn[n]
        end
        for n in indices
            ############### what it does ###############
            # Lr[m,n] = (s*sn[n])*(ϕ(abs(sum - Lrn[n])))
            ############################################
            if @inbounds sn[n] == s
                @inbounds Lr[m,n] = -ϕ(abs(sum - Lrn[n]))
            else
                @inbounds Lr[m,n] = ϕ(abs(sum - Lrn[n]))
            end
        end    
    end
end

# TABLE
function
    llr_horizontal_update!(
        Lr::Matrix{Float64},                            
        Lq::Matrix{Float64},
        indices_row::Vector{Vector{Int64}},
        sn::Vector{Bool},
        Lrn::Vector{Float64},
        phi::Vector{Float64}
    )

    m = 0
    for indices in indices_row
        m += 1
        sum = 0.0
        s = 1 
        for n in indices
            ####### what it does: #######
            # Lrn[n] = ϕ(abs(Lq[n,m]))
            # sum += Lrn[n]
            # sn[n] = sign(Lq[n,m])
            # s *= sn[n]
            #############################
            if @inbounds Lq[n,m] >= 0
                @inbounds Lrn[n] = phi[get_index(Lq[n,m])]
                @inbounds sn[n] = false
            else
                @inbounds Lrn[n] = phi[get_index(-Lq[n,m])]
                @inbounds sn[n] = true
                @inbounds s ⊻= sn[n]
            end
            @inbounds @fastmath sum += Lrn[n]
        end
        for n in indices
            ############### what it does ###############
            # Lr[m,n] = (s*sn[n])*(ϕ(abs(sum - Lrn[n])))
            ############################################
            if @inbounds sn[n] == s
                @inbounds @fastmath Lr[m,n] = -phi[get_index(abs(sum - Lrn[n]))]
            else
                @inbounds @fastmath Lr[m,n] = phi[get_index(abs(sum - Lrn[n]))]
            end
        end    
    end
end