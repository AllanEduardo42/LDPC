################################################################################
# Allan Eduardo Feitosa
# 26 set 2024
# Horizontal update of the LLR
# There are 3 different methods for SPA: tanh, alternative and table

include("min_sum.jl")

######################### SPA USING HYPERBOLIC TANGENT #########################
function
    _llr_horizontal_update!(
        Lr::AbstractVector{<:AbstractFloat},
        Lq::AbstractVector{<:AbstractFloat},
        nodes::Vector{<:Integer},
        Lrn::Vector{<:AbstractFloat},
        sn::Nothing,
        phi::Nothing        
    )
    pLr = 1.0
    for node in nodes
        @inbounds @fastmath Lrn[node] = tanh(0.5*Lq[node])
        @inbounds @fastmath pLr *= Lrn[node]
    end
    for node in nodes
        @inbounds @fastmath x = pLr/Lrn[node]
        if abs(x) < 1 # controls divergent values of Lr
            @inbounds @fastmath Lr[node] = 2*atanh(x)
        end
    end
end



###################### ALTERNATIVE TO HYPERBOLIC TANGENT #######################

function ϕ(Lq::AbstractFloat)    
    @fastmath log(1 + 2/(exp(Lq)-1))
end

function phi_sign!(Lq::AbstractFloat,s::Integer)
    return ϕ(abs(Lq)), sign(Lq), flipsign(s,Lq)
end

function 
    _llr_horizontal_update!(
        Lr::AbstractVector{<:AbstractFloat},
        Lq::AbstractVector{<:AbstractFloat},
        nodes::Vector{<:Integer},
        Lrn::Vector{<:AbstractFloat},
        sn::Vector{<:Integer},
        phi::Nothing
    )

    sLr = 0.0
    s = Int8(1) 
    for node in nodes
        @inbounds Lrn[node], sn[node], s = phi_sign!(Lq[check,node],s)
        @inbounds @fastmath sLr += Lrn[node] 
    end
    for node in nodes
        x = abs(sLr - Lrn[node])
        if x > 0 # (Inf restriction)
            @inbounds Lr[check,node] = flipsign(s,sn[node])*ϕ(x)
        end
    end    
end

############################ SPA USING LOOKUP TABLE ############################

function 
    phi_sign!(
        Lq::AbstractFloat,
        s::Integer,
        phi::Vector{<:AbstractFloat}
    )    
    return phi[get_index(abs(Lq))], sign(Lq), flipsign(s,Lq)
end

function
    _llr_horizontal_update!(
        Lr::AbstractVector{<:AbstractFloat},                            
        Lq::AbstractVector{<:AbstractFloat},
        nodes::Vector{<:Integer},
        Lrn::Vector{<:AbstractFloat},
        sn::Vector{Bool},
        phi::Vector{<:AbstractFloat}
    )
    
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

################################### MIN SUM ###################################


function
    _llr_horizontal_update!(
        Lr::AbstractVector{<:AbstractFloat},
        Lq::AbstractVector{<:AbstractFloat},
        nodes::Vector{<:Integer},
        Lrn::Nothing,
        sn::Vector{Bool},
        phi::Nothing        
    )

    args = _min_sum!(view(Lq,check,:),sn,nodes)
    for node in nodes
        Lr[check,node] = __min_sum!(node,sn[node],args...)
    end 

end

### The function below is used in the RBP algorithm
function
    llr_horizontal_update_one_check_only!(
        Lq::AbstractVector{<:AbstractFloat},
        nodes::Vector{<:Integer},
        _node::Integer,
        Lr::AbstractFloat
    )
    pLr = 1.0
    for node in nodes
        if node != _node
            @inbounds @fastmath pLr *= tanh(0.5*Lq[node])
        end
    end

    if abs(pLr) < 1 
        return @inbounds @fastmath 2*atanh(pLr)
    else
        return Lr
    end

end