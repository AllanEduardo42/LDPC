################################################################################
# Allan Eduardo Feitosa
# 22 set 2024
# Horizontal update of the LLR based MIN SUM Algorithm

function
    min_sum!(
        Lr::Matrix{<:AbstractFloat},                           
        Lq::Matrix{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        sn::Vector{Bool},
    )    

    check = 0
    for nodes in checks2nodes
        check += 1
        args = _min_sum!(
            view(Lq,check,:),
            sn,
            nodes
        )
        for node in nodes
            Lr[check,node] = __min_sum!(node,sn[node],args...)
        end
    end
end

function abs_sign!(Lq::AbstractFloat,s::Integer)
    sn = signbit(Lq)
    return abs(Lq), sn, s ⊻ sn
end

function _min_sum!(                           
    Lq::AbstractVector{<:AbstractFloat},
    sn::Vector{Bool},
    nodes::Vector{<:Integer},        
    )

    s = false
    minL = Inf
    minL2 = Inf
    max_idx = 0
    for node in nodes
        @inbounds β, sn[node], s = abs_sign!(Lq[node],s)
        if β < minL
            max_idx = node
            minL, minL2 = β, minL
        elseif β < minL2
            minL2 = β
        end
    end

    return minL, minL2, s, max_idx

end

function 
    __min_sum!(
        node::Integer,
        sn::Bool,
        minL::AbstractFloat,
        minL2::AbstractFloat,
        s::Bool,
        max_idx::Integer
    )

    if node == max_idx #(pick the second least Lq)
        if s ⊻ sn #if negative
            return -minL2
        else
            return minL2
        end
    else
        if s ⊻ sn #if negative
            return -minL
        else
            return minL
        end
    end

end

### specialized method for the RBP algorithm
function
    min_sum_RBP!(
        max_coords::Vector{<:Integer},
        max_residue::AbstractFloat,
        penalty::Matrix{<:AbstractFloat},
        Lr::AbstractVector{<:AbstractFloat},                           
        Lq::AbstractVector{<:AbstractFloat},
        sn::Vector{Bool},
        nodes::Vector{<:Integer},
        nmax::Integer,
        check::Integer
    )
    
    x = 0.0
    y = 0.0
    args = _min_sum!(Lq,sn,nodes)
    for node in nodes
        if node ≠ nmax
            x = __min_sum!(node,sn[node],args...)
            y = abs(x - Lr[node])*penalty[check,node]
            if y > max_residue
                max_residue = y
                max_coords[1] = check
                max_coords[2] = node
            end
        end
    end

    return max_residue
end

function
    min_sum_RBP_init!(          
        max_coords::Vector{<:Integer},              
        Lq::AbstractVector{<:AbstractFloat},
        sn::Vector{Bool},
        nodes::Vector{<:Integer},
        check::Integer
    )
    
    x = 0.0
    y = 0.0
    max_residue = 0.0
    args = _min_sum!(Lq,sn,nodes)
    for node in nodes
        x = __min_sum!(node,sn[node],args...)
        y = abs(x) 
        if y > max_residue
            max_residue = y
            max_coords[1] = check
            max_coords[2] = node
        end
    end
end

### specialized method for the RBP algorithm
function
    min_sum_RBP_R!(
        R::Matrix{<:AbstractFloat},
        Lr::AbstractVector{<:AbstractFloat},                           
        Lq::AbstractVector{<:AbstractFloat},
        sn::Vector{Bool},
        nodes::Vector{<:Integer},
        nmax::Integer,
        check::Integer
    )
    
    x = 0.0
    minL, minL2, s, max_idx = _min_sum!(Lq,sn,nodes)
    for node in nodes
        if node ≠ nmax
            x = __min_sum!(node,max_idx,s,sn[node],minL,minL2)
            R[check,node] = abs(x - Lr[node]) 
        end
    end

end

function
    min_sum_RBP_R_init!(     
        R::Matrix{<:AbstractFloat},                    
        Lq::AbstractVector{<:AbstractFloat},
        sn::Vector{<:Integer},
        nodes::Vector{<:Integer},
        check::Integer,
    )
    
    x = 0.0
    args = _min_sum!(Lq,sn,nodes)
    for node in nodes
        x = __min_sum!(node,sn[node],args...)
        R[check,node] = abs(x) 
    end

end

