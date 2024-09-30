################################################################################
# Allan Eduardo Feitosa
# 24 set 2024
# Specialized methods for the RBP algorithm

include("minsum.jl")

function
    minsum_RBP!(
        Residues::Nothing,
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
    args = _minsum!(Lq,signs,m,cn2vn)
    for n in cn2vn[m]
        if n ≠ vnmax
            x = __minsum!(n,signs[n],args...)
            y = abs(x - Lr[m,n])*Factors[m,n]
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
    minsum_RBP!(
        Residues::Matrix{<:AbstractFloat},
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
    args = _minsum!(Lq,signs,m,cn2vn)
    for n in cn2vn[m]
        if n ≠ vnmax
            x = __minsum!(n,signs[n],args...)
            Residues[m,n] = abs(x - Lr[m,n])*Factors[m,n]
        end
    end

    return maxresidue

end

function
    minsum_RBP_init!(    
        Residues::Nothing,      
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
        args = _minsum!(Lq,signs,m,cn2vn)
        for n in cn2vn[m]
            x = __minsum!(n,signs[n],args...)
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
    minsum_RBP_init!(     
        Residues::Matrix{<:AbstractFloat}, 
        maxcoords::Vector{<:Integer},                   
        Lq::Matrix{<:AbstractFloat},
        signs::Vector{<:Integer},
        cn2vn::Vector{Vector{T}} where {T<:Integer}
    )
    
    x = 0.0

    for m in eachindex(cn2vn)
        args = _minsum!(Lq,signs,m,cn2vn)
        for n in cn2vn[m]
            x = __minsum!(n,signs[n],args...)
            Residues[m,n] = abs(x) 
        end
    end
end