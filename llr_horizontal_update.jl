################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Horizontal update of the LLR based Sum-Product Algorithm

######################### SPA USING HYPERBOLIC TANGENT #########################
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
                if nn ≠ n
                    @inbounds @fastmath pLr *= tanh(0.5*Lq[nn,m])
                end
            end
            @inbounds @fastmath Lr[m,n] = 2*atanh(pLr)
        end
    end
end

###################### ALTERNATIVE TO HYPERBOLIC TANGENT #######################

function ϕ(x::Float64)::Float64
    @fastmath log(1 + 2/(exp(x)-1))
end

function phi_sign!(x::Float64,s::Int8)
    return ϕ(abs(x)), sign(x), flipsign(s,x)
end

function 
    llr_horizontal_update!(
        Lr::Matrix{Float64},                            
        Lq::Matrix{Float64},
        indices_row::Vector{Vector{Int64}},
        sn::Vector{Int8},
        Lrn::Vector{Float64},
        x::Nothing
    )

    m = 0
    for indices in indices_row
        m += 1
        sum = 0.0
        s = Int8(1) 
        for n in indices
            @inbounds Lrn[n], sn[n], s = phi_sign!(Lq[n,m],s)
            @inbounds @fastmath sum += Lrn[n] 
        end
        for n in indices
            @inbounds Lr[m,n] = flipsign(s,sn[n])*ϕ(abs(sum - Lrn[n]))
        end    
    end
end

############################ SPA USING LOOKUP TABLE ############################

function 
    phi_sign!(
        x::Float64,
        s::Int8,
        phi::Vector{Float64}
    )    
    return phi[get_index(abs(x))], sign(x), flipsign(s,x)
end

function
    llr_horizontal_update!(
        Lr::Matrix{Float64},                            
        Lq::Matrix{Float64},
        indices_row::Vector{Vector{Int64}},
        sn::Vector{Int8},
        Lrn::Vector{Float64},
        phi::Vector{Float64}
    )

    m = 0
    for indices in indices_row
        m += 1
        sum = 0.0
        s = Int8(1) 
        for n in indices
            @inbounds Lrn[n], sn[n], s = phi_sign!(Lq[n,m],s,phi)
            @inbounds @fastmath sum += Lrn[n]
        end
        for n in indices
            @inbounds Lr[m,n] = flipsign(s,sn[n])*phi[get_index(abs(sum - Lrn[n]))]
        end    
    end
end

################################### MIN SUM ####################################

function abs_sign!(x::Float64,s::Int8)
    return abs(x), sign(x), flipsign(s,x)
end

function
    llr_horizontal_update!(
        Lr::Matrix{Float64},                           
        Lq::Matrix{Float64},
        indices_row::Vector{Vector{Int64}},
        sn::Vector{Int8},
        x::Nothing,
        y::Nothing
    )    

    m = 0
    for indices in indices_row
        m += 1       
        s = Int8(1) 
        # first index
        @inbounds idx = indices[1]
        @inbounds minL, sn[idx], s = abs_sign!(Lq[idx,m],s) 
        # second index
        @inbounds n = indices[2]
        @inbounds β, sn[n], s = abs_sign!(Lq[n,m],s)      
        if β < minL
            idx = n
            minL2 = minL
            minL = β
        else
            minL2 = β
        end
        # remaining indices
        next = iterate(indices,3)
        while next !== nothing
            (n,state) = next
            @inbounds β, sn[n], s = abs_sign!(Lq[n,m],s)
            if β < minL
                idx = n
                minL2 = minL
                minL = β
            elseif β < minL2
                minL2 = β
            end
            next = iterate(indices,state)
        end
        for n in indices
            if n == idx
                @inbounds Lr[m,n] = flipsign(s,sn[n])*minL2
            else
                @inbounds Lr[m,n] = flipsign(s,sn[n])*minL
            end
        end  
    end
end