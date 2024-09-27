################################################################################
# Allan Eduardo Feitosa
# 24 set 2024
# Specialized methods for the RBP algorithm

include("min_sum.jl")

function
    min_sum_lRBP!(
        max_coords::Vector{<:Integer},
        max_residue::AbstractFloat,
        penalty::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},                           
        Lq::Matrix{<:AbstractFloat},
        sn::Vector{Bool},
        nodes::Vector{<:Integer},
        nmax::Integer,
        check::Integer
    )
    
    x = 0.0
    y = 0.0
    args = _min_sum!(Lq,check,sn,nodes)
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
    min_sum_lRBP_init!(          
        max_coords::Vector{<:Integer},              
        Lq::Matrix{<:AbstractFloat},
        sn::Vector{Bool},
        checks2nodes::Vector{Vector{T}} where {T<:Integer}
    )
    
    x = 0.0
    y = 0.0
    max_residue = 0.0
    check = 0
    for nodes in checks2nodes
        check += 1
        args = _min_sum!(Lq,check,sn,nodes)
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
end

### specialized method for the RBP algorithm
function
    min_sum_RBP!(
        Residues::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},                           
        Lq::Matrix{<:AbstractFloat},
        sn::Vector{Bool},
        nodes::Vector{<:Integer},
        nmax::Integer,
        check::Integer
    )
    
    x = 0.0
    args = _min_sum!(Lq,check,sn,nodes)
    for node in nodes
        if node ≠ nmax
            x = __min_sum!(node,sn[node],args...)
            Residues[check,node] = abs(x - Lr[check,node]) 
        end
    end

end

function
    min_sum_RBP_init!(     
        Residues::Matrix{<:AbstractFloat},                    
        Lq::Matrix{<:AbstractFloat},
        sn::Vector{<:Integer},
        checks2nodes::Vector{Vector{T}} where {T<:Integer}
    )
    
    x = 0.0
    check = 0
    for nodes in checks2nodes
        check += 1
        args = _min_sum!(Lq,check,sn,nodes)
        for node in nodes
            x = __min_sum!(node,sn[node],args...)
            Residues[check,node] = abs(x) 
        end
    end
end