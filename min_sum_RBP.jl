################################################################################
# Allan Eduardo Feitosa
# 24 set 2024
# Specialized methods for the RBP algorithm

include("min_sum.jl")

function
    min_sum_lRBP!(
        maxcoords::Vector{<:Integer},
        maxresidue::AbstractFloat,
        Factors::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},                           
        Lq::Matrix{<:AbstractFloat},
        signs::Vector{Bool},
        vnmax::Integer,
        m::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer}       
    )
    
    x = 0.0
    y = 0.0
    args = _min_sum!(Lq,m,signs,cn2vn[m])
    for n in cn2vn[m]
        if n ≠ vnmax
            x = __min_sum!(n,signs[n],args...)
            y = abs(x - Lr[n])*Factors[m,n]
            if y > maxresidue
                maxresidue = y
                maxcoords[1] = m
                maxcoords[2] = n
            end
        end
    end

    return maxresidue
end

function
    min_sum_lRBP_init!(          
        maxcoords::Vector{<:Integer},              
        Lq::Matrix{<:AbstractFloat},
        signs::Vector{Bool},
        cn2vn::Vector{Vector{T}} where {T<:Integer}
    )
    
    x = 0.0
    y = 0.0
    maxresidue = 0.0
    for m in eachindex(cn2vn)
        cn2vn[m] = cn2vn[m]
        args = _min_sum!(Lq,signs,m,cn2vn)
        for n in cn2vn[m]
            x = __min_sum!(n,signs[n],args...)
            y = abs(x) 
            if y > maxresidue
                maxresidue = y
                maxcoords[1] = m
                maxcoords[2] = n
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
        signs::Vector{Bool},
        vnmax::Integer,
        m::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer}
    )
    
    x = 0.0
    args = _min_sum!(Lq,signs,m,cn2vn)
    for n in cn2vn[m]
        if n ≠ vnmax
            x = __min_sum!(n,signs[n],args...)
            Residues[m,n] = abs(x - Lr[m,n]) 
        end
    end

end

function
    min_sum_RBP_init!(     
        Residues::Matrix{<:AbstractFloat},                    
        Lq::Matrix{<:AbstractFloat},
        signs::Vector{<:Integer},
        cn2vn::Vector{Vector{T}} where {T<:Integer}
    )
    
    x = 0.0

    for m in eachindex(cn2vn)
        args = _min_sum!(Lq,signs,m,cn2vn)
        for n in cn2vn[m]
            x = __min_sum!(n,signs[n],args...)
            Residues[m,n] = abs(x) 
        end
    end
end