################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Horizontal update of the LLR based Sum-Product Algorithm (with Inf restriction)

include("tanh_llr_horizontal_update.jl")

######################### SPA USING HYPERBOLIC TANGENT #########################
function 
    llr_horizontal_update!(
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        Lrn::Vector{<:AbstractFloat},
        x::Nothing,
        z::Nothing
    )

    check = 0
    for nodes in checks2nodes
        check += 1
        tanh_llr_horizontal_update!(Lr,Lq,nodes,check,Lrn)
    end
end

###################### ALTERNATIVE TO HYPERBOLIC TANGENT #######################

function ϕ(x::AbstractFloat)    
    @fastmath log(1 + 2/(exp(x)-1))
end

function phi_sign!(x::AbstractFloat,s::Integer)
    return ϕ(abs(x)), sign(x), flipsign(s,x)
end

function 
    llr_horizontal_update!(
        Lr::Matrix{<:AbstractFloat},                            
        Lq::Matrix{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        Lrn::Vector{<:AbstractFloat},
        sn::Vector{<:Integer},
        x::Nothing
    )

    check = 0
    for nodes in checks2nodes
        check += 1
        sLr = 0.0
        s = 0x01 
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
        x::AbstractFloat,
        s::Integer,
        phi::Vector{<:AbstractFloat}
    )    
    return phi[get_index(abs(x))], sign(x), flipsign(s,x)
end

function
    llr_horizontal_update!(
        Lr::Matrix{<:AbstractFloat},                            
        Lq::Matrix{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        Lrn::Vector{<:AbstractFloat},
        sn::Vector{<:Integer},
        phi::Vector{<:AbstractFloat}
    )

    check = 0
    for nodes in checks2nodes
        check += 1
        sLr = 0.0
        s = 0x01
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

function abs_sign!(x::AbstractFloat,s::Integer)
    return abs(x), sign(x), flipsign(s,x)
end

function
    llr_horizontal_update!(
        Lr::Matrix{<:AbstractFloat},                           
        Lq::Matrix{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        x::Nothing,
        sn::Vector{<:Integer},
        y::Nothing
    )    

    check = 0
    for nodes in checks2nodes
        check += 1       
        s = Integer(1) 
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