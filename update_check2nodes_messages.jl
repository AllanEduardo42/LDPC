################################################################################
# Allan Eduardo Feitosa
# 26 set 2024
# Horizontal update of the LLR
# There are 4 different methods for SPA: Mckay, tanh, alternative and table

include("min_sum.jl")

########################### SPA USING MKAY's METHOD ############################
function
    update_check2nodes_messages!(
        r::Array{<:AbstractFloat,3},
        δq::Matrix{<:AbstractFloat},
        check::Integer,
        nodes::Vector{<:Integer},     
    )
    δr = 1
    for node in nodes
        @inbounds δr *= δq[check,node]
    end
    for node in nodes
        x = δr/(δq[check,node] + eps())        
        @inbounds r[check,node,1] = 0.5*(1+x)
        @inbounds r[check,node,2] = 0.5*(1-x)
    end
end

######################### SPA USING HYPERBOLIC TANGENT #########################
function
    update_check2nodes_messages!(
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        check::Integer,
        nodes::Vector{<:Integer},
        Lrn::Vector{<:AbstractFloat},
        sn::Nothing,
        phi::Nothing        
    )
    pLr = 1.0
    for node in nodes
        @inbounds @fastmath Lrn[node] = tanh(0.5*Lq[check,node])
        @inbounds @fastmath pLr *= Lrn[node]
    end
    for node in nodes
        @inbounds @fastmath x = pLr/Lrn[node]
        if abs(x) < 1 # controls divergent values of Lr
            @inbounds @fastmath Lr[check,node] = 2*atanh(x)
        end
    end
end



################# ALTERNATIVE TO HYPERBOLIC TANGENT AND TABLE ##################

function ϕ(Lq::AbstractFloat, phi::Nothing)    
    @fastmath log(1 + 2/(exp(Lq)-1))
end

function ϕ(Lq::AbstractFloat, phi::Vector{<:AbstractFloat})    
    phi[get_index(Lq)]
end

function
    phi_sign!(
        Lq::AbstractFloat,
        s::Bool,
        phi::Union{Vector{<:AbstractFloat},Nothing}
    )

    sn = signbit(Lq)
    return ϕ(abs(Lq),phi), sn, s ⊻ sn

end

function 
    update_check2nodes_messages!(
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        check::Integer,
        nodes::Vector{<:Integer},
        Lrn::Vector{<:AbstractFloat},
        sn::Vector{Bool},
        phi::Union{Vector{<:AbstractFloat},Nothing}
    )

    sLr = 0.0
    s = false
    for node in nodes
        @inbounds Lrn[node], sn[node], s = phi_sign!(Lq[check,node],s,phi)
        @inbounds @fastmath sLr += Lrn[node] 
    end
    for node in nodes
        x = abs(sLr - Lrn[node])
        if x > 0 # (Inf restriction)
            @inbounds Lr[check,node] = (1 - 2*(sn[node] ⊻ s))*ϕ(x,phi)
        end
    end    
end

################################### MIN SUM ###################################


function
    update_check2nodes_messages!(
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        check::Integer,
        nodes::Vector{<:Integer},
        Lrn::Nothing,
        sn::Vector{Bool},
        phi::Nothing        
    )
    args = _min_sum!(Lq,check,sn,nodes)
    for node in nodes
        Lr[check,node] = __min_sum!(node,sn[node],args...)
    end 

end

### The function below is used in the RBP algorithm
function
    update_check2node_message!(
        Lq::Matrix{<:AbstractFloat},
        check::Integer,
        nodes::Vector{<:Integer},
        _node::Integer,
        Lr::AbstractFloat
    )
    pLr = 1.0
    for node in nodes
        if node != _node
            @inbounds @fastmath pLr *= tanh(0.5*Lq[check,node])
        end
    end

    if abs(pLr) < 1 
        return @inbounds @fastmath 2*atanh(pLr)
    else
        return Lr
    end

end