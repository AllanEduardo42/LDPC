################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Horizontal update of the LLR based Sum-Product Algorithm

######################### SPA USING HYPERBOLIC TANGENT #########################
### No Inf restriction ###
function 
    llr_horizontal_update!(
        Lr::Matrix{Float64},
        Lq::Matrix{Float64},
        checks2nodes::Vector{Vector{Int64}},
        x::Nothing,
        y::Nothing,
        z::Nothing
    )

    check = 0
    for nodes in checks2nodes
        check += 1
        for node in nodes
            @inbounds Lr[check,node] = 1.0
            for n in nodes
                if n ≠ node
                    @inbounds @fastmath Lr[check,node] *= tanh(0.5*Lq[check,n])
                end
            end
            @inbounds @fastmath Lr[check,node] = 2*atanh(Lr[check,node])
        end
    end
end

### Inf restriction ###
function 
    llr_horizontal_update!(
        Lr::Matrix{Float64},
        Lq::Matrix{Float64},
        checks2nodes::Vector{Vector{Int64}},
        x::Nothing,
        Lrn::Vector{Float64},
        z::Nothing
    )

    check = 0
    for nodes in checks2nodes
        check += 1
        pLr = 1.0
        for node in nodes
            @inbounds @fastmath Lrn[node] = tanh(0.5*Lq[check,node])
            @inbounds @fastmath pLr *= Lrn[node]
        end
        for node in nodes
            @inbounds @fastmath x = pLr/Lrn[node]
            if abs(x) < 1
                @inbounds @fastmath Lr[check,node] = 2*atanh(x)
            end
        end
    end
end

###################### ALTERNATIVE TO HYPERBOLIC TANGENT #######################

function ϕ(x::Float64)    
    @fastmath log(1 + 2/(exp(x)-1))
end

function phi_sign!(x::Float64,s::Int8)
    return ϕ(abs(x)), sign(x), flipsign(s,x)
end

function 
    llr_horizontal_update!(
        Lr::Matrix{Float64},                            
        Lq::Matrix{Float64},
        checks2nodes::Vector{Vector{Int64}},
        sn::Vector{Int8},
        Lrn::Vector{Float64},
        x::Nothing
    )

    check = 0
    for nodes in checks2nodes
        check += 1
        sLr = 0.0
        s = Int8(1) 
        for node in nodes
            @inbounds Lrn[node], sn[node], s = phi_sign!(Lq[check,node],s)
            @inbounds @fastmath sLr += Lrn[node] 
        end
        for node in nodes
            x = abs(sLr - Lrn[node])
            if x > 0 
                @inbounds Lr[check,node] = flipsign(s,sn[node])*ϕ(x)
            end
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
        checks2nodes::Vector{Vector{Int64}},
        sn::Vector{Int8},
        Lrn::Vector{Float64},
        phi::Vector{Float64}
    )

    check = 0
    for nodes in checks2nodes
        check += 1
        sLr = 0.0
        s = Int8(1) 
        for node in nodes
            @inbounds Lrn[node], sn[node], s = phi_sign!(Lq[check,node],s,phi)
            @inbounds @fastmath sLr += Lrn[node]
        end
        for node in nodes
            @inbounds Lr[check,node] = flipsign(s,sn[node])*phi[get_index(sLr - Lrn[node])]
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
        checks2nodes::Vector{Vector{Int64}},
        sn::Vector{Int8},
        x::Nothing,
        y::Nothing
    )    

    check = 0
    for nodes in checks2nodes
        check += 1       
        s = Int8(1) 
        # first index
        @inbounds idx = nodes[1]
        @inbounds minL, sn[idx], s = abs_sign!(Lq[check,idx],s) 
        # second index
        @inbounds node = nodes[2]
        @inbounds β, sn[node], s = abs_sign!(Lq[check,node],s)      
        if β < minL
            idx = node
            minL2 = minL
            minL = β
        else
            minL2 = β
        end
        # remaining nodes
        next = iterate(nodes,3)
        while next !== nothing
            (node,state) = next
            @inbounds β, sn[node], s = abs_sign!(Lq[check,node],s)
            if β < minL
                idx = node
                minL2 = minL
                minL = β
            elseif β < minL2
                minL2 = β
            end
            next = iterate(nodes,state)
        end
        for node in nodes
            if node == idx
                @inbounds Lr[check,node] = flipsign(s,sn[node])*minL2
            else
                @inbounds Lr[check,node] = flipsign(s,sn[node])*minL
            end
        end  
    end
end