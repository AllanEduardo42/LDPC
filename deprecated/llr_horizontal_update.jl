################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Horizontal update of the LLR based Sum-Product Algorithm
# There are 3 different methods to update: tanh, alternative and table

######################### SPA USING HYPERBOLIC TANGENT #########################
function 
    llr_horizontal_update_tnh!(
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        Lrn::Vector{<:AbstractFloat}
    )

    check = 0
    for nodes in checks2nodes
        check += 1
        _llr_horizontal_update_tnh!(
            view(Lr,check,:),
            view(Lq,check,:),
            Lrn,
            nodes,
        )
    end
end

# The core function below is is also used in the LBP algorithm
function
    _llr_horizontal_update_tnh!(
        Lr::AbstractVector{<:AbstractFloat},
        Lq::AbstractVector{<:AbstractFloat},
        Lrn::Vector{<:AbstractFloat},
        nodes::Vector{<:Integer},        
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

###################### ALTERNATIVE TO HYPERBOLIC TANGENT #######################

function ϕ(Lq::AbstractFloat)    
    @fastmath log(1 + 2/(exp(Lq)-1))
end

function phi_sign!(Lq::AbstractFloat,s::Integer)
    return ϕ(abs(Lq)), sign(Lq), flipsign(s,Lq)
end

function 
    llr_horizontal_update_alt!(
        Lr::Matrix{<:AbstractFloat},                            
        Lq::Matrix{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        Lrn::Vector{<:AbstractFloat},
        sn::Vector{<:Integer},
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
            if x > 0 # (Inf restriction)
                @inbounds Lr[check,node] = flipsign(s,sn[node])*ϕ(x)
            end
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
    llr_horizontal_update_tab!(
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

