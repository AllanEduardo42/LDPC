################################################################################
# Allan Eduardo Feitosa
# 22 set 2024
# Horizontal update of the LLR based MIN SUM Algorithm

function abs_sign!(Lq::AbstractFloat,s::Integer)
    return abs(Lq), sign(Lq), flipsign(s,Lq)
end

function
    min_sum!(
        Lr::Matrix{<:AbstractFloat},                           
        Lq::Matrix{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        sn::Vector{<:Integer},
    )    

    check = 0
    for nodes in checks2nodes
        check += 1
        minL, minL2, s, max_idx = _min_sum!(view(Lq,check,:),sn,nodes)
        for node in nodes
            Lr[check,node] = __min_sum!(node,max_idx,s,sn[node],minL,minL2)
        end
    end
end

function _min_sum!(                           
    Lq::AbstractVector{<:AbstractFloat},
    sn::Vector{<:Integer},
    nodes::Vector{<:Integer},        
    )

    s = Int8(1)
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
        max_idx::Integer,
        s::Int8,
        sn::Int8,
        minL::AbstractFloat,
        minL2::AbstractFloat
    )

    if node == max_idx #(pick the second least Lq)
        return flipsign(s,sn)*minL2
    else
        return flipsign(s,sn)*minL
    end

end

### specialized method for the RBP algorithm
function
    min_sum_RBP!(
        Lr::AbstractVector{<:AbstractFloat},                           
        Lq::AbstractVector{<:AbstractFloat},
        sn::Vector{<:Integer},
        nodes::Vector{<:Integer},
        nmax::Integer,
        check::Integer,
        max_residue::AbstractFloat,
        max_coords::Vector{<:Integer} 
    )
    
    x = 0.0
    y = 0.0
    minL, minL2, s, max_idx = _min_sum!(Lq,sn,nodes)
    for node in nodes
        if node ≠ nmax
            x = __min_sum!(node,max_idx,s,sn[node],minL,minL2)
            y = abs(x - Lr[node]) 
            if y > max_residue
                max_residue = y
                max_coords[1] = check
                max_coords[2] = node
            end
        end
    end

    return max_residue
end

### specialized method for the RBP algorithm
function
    min_sum_RBP_R!(
        Lr::AbstractVector{<:AbstractFloat},                           
        Lq::AbstractVector{<:AbstractFloat},
        sn::Vector{<:Integer},
        nodes::Vector{<:Integer},
        nmax::Integer,
        check::Integer,
        R::Matrix{<:AbstractFloat}
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
    min_sum_RBP_R!(
        Lr::AbstractVector{<:AbstractFloat},                           
        Lq::AbstractVector{<:AbstractFloat},
        sn::Vector{<:Integer},
        nodes::Vector{<:Integer},
        check::Integer,
        R::Matrix{<:AbstractFloat}
    )
    
    x = 0.0
    minL, minL2, s, max_idx = _min_sum!(Lq,sn,nodes)
    for node in nodes
        x = __min_sum!(node,max_idx,s,sn[node],minL,minL2)
        R[check,node] = abs(x - Lr[node]) 
    end

end

