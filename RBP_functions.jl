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
        listres::Union{Vector{<:AbstractFloat},Nothing},
        listadd::Union{Matrix{<:Integer},Nothing},
        listsize::Integer,
        inlist::Union{Matrix{<:Integer},Nothing}           
    )
    
    minsum!(Lq,Ms,signs,m,cn2vn)
    x = 0.0    
    @inbounds for n in cn2vn[m]
        if n ≠ vnmax
            x = calc_residue(Ms,Factors,Lr,m,n)
            maxresidue = findmaxresidue!(Residues,maxcoords,maxresidue,m,n,x,
                listres,listadd,listsize,inlist)
            
        end
    end
    if iszero(maxcoords) || maxresidue == 0.0
        maxcoords[1] = rand(1:length(cn2vn))
        maxcoords[2] = rand(cn2vn[maxcoords[1]])
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
        maxcoords::Vector{<:Integer},
        maxresidue::AbstractFloat,
        m::Integer,
        n::Integer,
        x::AbstractFloat,
        ::Nothing,
        ::Nothing,
        ::Integer,
        ::Nothing
    )

    @inbounds Residues[m,n] = x

    if @fastmath x > maxresidue
        maxresidue = x
        @inbounds maxcoords[1] = m
        @inbounds maxcoords[2] = n
    end

    return maxresidue
end

# Local-RBP
function 
    findmaxresidue!(
        ::Nothing,
        maxcoords::Vector{<:Integer},
        maxresidue::AbstractFloat,
        m::Integer,
        n::Integer,
        x::AbstractFloat,
        ::Nothing,
        ::Nothing,
        ::Integer,
        ::Nothing
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
        maxcoords::Vector{<:Integer},
        maxresidue::AbstractFloat,
        m::Integer,
        n::Integer,
        x::AbstractFloat,
        listres::Vector{<:AbstractFloat},
        listadd::Matrix{<:Integer},
        listsize::Integer,
        inlist::Matrix{<:Integer}
    )

    if @inbounds inlist[m,n]
        @inbounds inlist[m,n] = false
        @inbounds for i in 1:listsize
            mm = listadd[1,i]
            nn = listadd[2,i]
            if mm == m && nn == n                
                for j=i:listsize
                    listres[j] = listres[j+1]
                    listadd[1,j] = listadd[1,j+1]
                    listadd[2,j] = listadd[2,j+1]
                end
                break
            end
        end
    end

    @fastmath @inbounds for i in 1:listsize
        y = listres[i]
        if x > y
            mm = listadd[1,end]
            nn = listadd[2,end]
            if mm ≠ 0
                inlist[nn,nn] = false
            end
            for j=listsize:-1:i+1
                listres[j] = listres[j-1]
                listadd[1,j] = listadd[1,j-1]
                listadd[2,j] = listadd[2,j-1]
            end
            listadd[1,i] = m
            listadd[2,i] = n
            listres[i] = x
            inlist[m,n] = true
            break
        end
    end

    @inbounds maxcoords[1] = listadd[1,1]
    @inbounds maxcoords[2] = listadd[2,1]

    @inbounds return listres[1]
end

function
    init_residues!(    
        Residues::Union{Matrix{<:AbstractFloat},Nothing},      
        maxcoords::Vector{<:Integer},              
        Lq::Matrix{<:AbstractFloat},
        signs::Vector{Bool},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        Ms::Matrix{<:AbstractFloat},
        listres::Union{Vector{<:AbstractFloat},Nothing},
        listadd::Union{Matrix{<:Integer},Nothing},
        listsize::Integer,
        inlist::Union{Matrix{<:Integer},Nothing}         
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
            listres,
            listadd,
            listsize,
            inlist
        )
    end

    return maxresidue
end

# RBP
function
    find_maxresidue_coords!(
        maxresidue::AbstractFloat,
        maxcoords::Vector{<:Integer},
        Residues::Matrix{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        ::Nothing,
        ::Nothing
    )

    @fastmath @inbounds for m in eachindex(cn2vn)
        for n in cn2vn[m]
            residue = Residues[m,n]
            if residue > maxresidue
                maxresidue = residue
                maxcoords[1] = m
                maxcoords[2] = n
            end
        end
    end

    return maxresidue
end

# Random-RBP
function
    find_maxresidue_coords!(
        maxresidue::AbstractFloat,
        maxcoords::Vector{<:Integer},
        Residues::Matrix{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        samples::Vector{<:Integer},
        rng_sample::AbstractRNG
    )

    M = length(cn2vn)
    rand!(rng_sample,samples,1:M)
    @fastmath @inbounds for m in samples
        for n in cn2vn[m]
            residue = Residues[m,n]
            if residue > maxresidue
                maxresidue = residue
                maxcoords[1] = m
                maxcoords[2] = n
            end
        end
    end

    return maxresidue
end

# Local-RBP and List-RBP
function 
    find_maxresidue_coords!(
        maxresidue::AbstractFloat,
        ::Vector{<:Integer},
        ::Nothing,
        ::Vector{Vector{T}} where {T<:Integer},
        ::Nothing,
        ::Nothing
    )

    return maxresidue
end

