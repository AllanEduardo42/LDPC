################################################################################
# Allan Eduardo Feitosa
# 24 set 2024
# Specialized methods for the RBP algorithm

include("minsum.jl")

function
    calc_residues!(
        Residues::Union{Matrix{<:AbstractFloat},Nothing},
        maxcoords::Vector{<:Integer},
        maxresidue::AbstractFloat,
        Factors::Union{Matrix{<:AbstractFloat},Nothing},
        Ms::Matrix{<:AbstractFloat},
        Lr::Union{Matrix{<:AbstractFloat},Nothing},                      
        Lq::Matrix{<:AbstractFloat},
        signs::Vector{Bool},
        vnmax::Integer,
        m::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        list::Union{Vector{Tuple{Float64,Vector{Int}}},Nothing},
        listsize::Integer             
    )
    
    minsum!(Lq,Ms,signs,m,cn2vn)
    x = 0.0    
    for n in cn2vn[m]
        if n ≠ vnmax
            x = calc_residue(Ms,Factors,Lr,m,n)
            maxresidue = findmaxresidue!(Residues,maxcoords,maxresidue,m,n,x,
                list,listsize)
            
        end
    end
    return maxresidue
end

function 
    calc_residue(
        Ms::Matrix{<:AbstractFloat},
        Factors::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        m::Integer,
        n::Integer
    )

    @fastmath @inbounds return abs(Ms[m,n] - Lr[m,n])*Factors[m,n]

end

# for initialization (Lr[m,n] = 0.0 and Factors[m,n] = 1.0 ∀m,n)
function 
    calc_residue(
        Ms::Matrix{<:AbstractFloat},
        ::Nothing,
        ::Nothing,
        m::Integer,
        n::Integer
    )

    @fastmath @inbounds return abs(Ms[m,n])

end

# RBP and Random-RBP
function 
    findmaxresidue!(
        Residues::Matrix{<:AbstractFloat},
        ::Vector{Int},
        maxresidue::AbstractFloat,
        m::Integer,
        n::Integer,
        x::AbstractFloat,
        ::Nothing,
        ::Integer
    )

    @inbounds Residues[m,n] = x

    return maxresidue
end
# LRBP
function 
    findmaxresidue!(
        ::Nothing,
        maxcoords::Vector{Int},
        maxresidue::AbstractFloat,
        m::Integer,
        n::Integer,
        x::AbstractFloat,
        ::Nothing,
        ::Integer
    )
    if @fastmath x > maxresidue
        maxresidue = x
        @inbounds maxcoords[1] = m
        @inbounds maxcoords[2] = n
    end

    return maxresidue
end
# List-RBP
function 
    findmaxresidue!(
        ::Nothing,
        coords::Vector{Int},
        maxresidue::AbstractFloat,
        m::Integer,
        n::Integer,
        x::AbstractFloat,
        list::Vector{Tuple{Float64,Vector{Int}}},
        listsize::Integer
    )

    for i in 1:listsize
        y = @inbounds list[i][1]
        if @fastmath x > y
            for j=listsize:-1:i+1
                @inbounds list[j] = list[j-1]
            end
            @inbounds coords[1] = m
            @inbounds coords[2] = n
            @inbounds list[i] = (x,coords)
            break
        end
    end

    return maxresidue
end

function
    init_residues!(    
        Residues::Union{Matrix{<:AbstractFloat},Nothing},      
        maxcoords::Vector{<:Integer},              
        Lq::Matrix{<:AbstractFloat},
        signs::Vector{Bool},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        Ms::Matrix{<:AbstractFloat}        
    )
    
    maxresidue = 0.0
    for m in eachindex(cn2vn)
        maxresidue = calc_residues!(
            Residues,
            maxcoords,
            maxresidue,
            nothing,
            Ms,
            nothing,
            Lq,
            signs,
            0,
            m,
            cn2vn,
            nothing,
            0,
        )
    end
end


function
    find_maxresidue_coords!(
        maxresidue::AbstractFloat,
        maxcoords::Vector{<:Integer},
        Residues::Matrix{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        samples::Vector{<:Integer},
        rng_sample::AbstractRNG,
        list::Nothing,
        listsize::Nothing
    )

    M = length(cn2vn)
    rand!(rng_sample,samples,1:M)
    for m in samples
        for n in cn2vn[m]
            @inbounds residue = Residues[m,n]
            if @fastmath residue > maxresidue
                maxresidue = residue
                @inbounds maxcoords[1] = m
                @inbounds maxcoords[2] = n
            end
        end
    end

    return maxresidue
end

function
    find_maxresidue_coords!(
        maxresidue::AbstractFloat,
        maxcoords::Vector{<:Integer},
        Residues::Matrix{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        ::Nothing,
        ::Nothing,
        list::Nothing,
        listsize::Nothing,
    )

    for m in eachindex(cn2vn)
        for n in cn2vn[m]
            @inbounds residue = Residues[m,n]
            if @fastmath residue > maxresidue
                maxresidue = residue
                @inbounds maxcoords[1] = m
                @inbounds maxcoords[2] = n
            end
        end
    end

    return maxresidue
end

