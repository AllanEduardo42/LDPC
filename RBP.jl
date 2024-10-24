################################################################################
# Allan Eduardo Feitosa
# 11 out 2024
# RBP Sum-Product Algorithm using min-sum to calculate the residues

include("RBP_functions.jl")

#RBP
function
    RBP!(
        Residues::Union{Matrix{<:AbstractFloat},Nothing},
        bitvector::Vector{Bool},
        Lr::Matrix{<:AbstractFloat},
        Ms::Matrix{<:AbstractFloat},
        maxresidue::AbstractFloat,
        maxcoords::Vector{<:Integer},
        Lq::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        signs::Vector{Bool},
        Factors::Matrix{<:AbstractFloat},
        rbpfactor::AbstractFloat,
        num_edges::Integer,
        Ldn::Vector{<:AbstractFloat},
        samples::Union{Vector{<:Integer},Nothing},
        rng_sample::Union{AbstractRNG,Nothing},
        listsize::Integer,
        listres::Union{Vector{<:AbstractFloat},Nothing},
        listadd::Union{Matrix{<:Integer},Nothing},
        inlist::Union{Matrix{<:Integer},Nothing}
    )

    for e in 1:num_edges

        if maxresidue == 0.0 # if RBP has converged
            break
        end      

        _RBP_update_Lr!(maxcoords,Factors,rbpfactor,cn2vn,Lq,Lr)

        # @inbounds Residues[maxcoords[1],maxcoords[2]] = 0.0
        __RBP(Residues,maxcoords,listsize,listres,listadd,inlist)

        maxresidue = _RBP_update_vn2cn!(Residues,maxcoords,0.0,Factors,Lf,Ldn,
            bitvector,vn2cn,cn2vn,Ms,Lr,Lq,signs,listsize,listres,listadd,inlist)

        maxresidue = find_maxresidue_coords!(maxresidue,maxcoords,Residues,cn2vn,
            samples,rng_sample)

    end
end

# RBP and Random-RBP
function 
    __RBP(
        Residues::Matrix{<:AbstractFloat},
        maxcoords::Vector{<:Integer},
        ::Integer,
        ::Nothing,
        ::Nothing,
        ::Nothing
    )

    @inbounds Residues[maxcoords[1],maxcoords[2]] = 0.0
end

# Local-RBP
function 
    __RBP(
        ::Nothing,
        ::Vector{<:Integer},
        ::Integer,
        ::Nothing,
        ::Nothing,
        ::Nothing
    )
    
end

# List-RBP
function 
    __RBP(
        ::Nothing,
        maxcoords::Vector{<:Integer},
        listsize::Integer,
        listres::Vector{<:AbstractFloat},
        listadd::Matrix{<:Integer},
        inlist::Matrix{Bool}
    )

    @inbounds inlist[maxcoords[1],maxcoords[2]] = false
    for i in 1:listsize   
        @inbounds listres[i] = listres[i+1]
        @inbounds listadd[1,i] = listadd[1,i+1]
        @inbounds listadd[2,i] = listadd[2,i+1]
    end
    
end

function 
    _RBP_update_Lr!(
        maxcoords::Vector{<:Integer},
        Factors::Matrix{<:AbstractFloat},
        rbpfactor::AbstractFloat,
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat}
    )

    (cnmax,vnmax) = maxcoords
    @fastmath @inbounds Factors[cnmax,vnmax] *= rbpfactor

    ### update Lr[cnmax,vnmax]
    pLr = 1.0
    for n in cn2vn[cnmax]
        if n != vnmax
            @fastmath @inbounds pLr *= tanh(0.5*Lq[n,cnmax])
        end
    end    
    if @fastmath abs(pLr) < 1 
        @fastmath @inbounds Lr[cnmax,vnmax] = 2*atanh(pLr)
    end

end

function 
    _RBP_update_vn2cn!(
        Residues::Union{Matrix{<:AbstractFloat},Nothing},
        maxcoords::Vector{<:Integer},
        maxresidue::AbstractFloat,
        Factors::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        Ldn::Vector{<:AbstractFloat},
        bitvector::Vector{Bool},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        Ms::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        signs::Vector{Bool},
        listsize::Integer,
        listres::Union{Vector{<:AbstractFloat},Nothing},
        listadd::Union{Matrix{<:Integer},Nothing},        
        inlist::Union{Matrix{<:Integer},Nothing}        
    )

    # update Ldn[vmax] and bitvector[vnmax]
    (cnmax,vnmax) = maxcoords
    @inbounds Ldn[vnmax] = Lf[vnmax]
    for m in vn2cn[vnmax]
        @fastmath @inbounds Ldn[vnmax] += Lr[m,vnmax]
        @fastmath @inbounds bitvector[vnmax] = signbit(Ldn[vnmax])
    end

    for m in vn2cn[vnmax]
        if m ≠ cnmax
            # update vn2cn messages Lq[vnmax,m], ∀m ≠ cnmax
            @fastmath @inbounds Lq[vnmax,m] = Ldn[vnmax] - Lr[m,vnmax]
            maxresidue = calc_residues!(Residues,maxcoords,maxresidue,Factors,
                                        Ms,Lr,Lq,signs,vnmax,m,cn2vn,listres,
                                        listadd,listsize,inlist)
        end
    end

    return maxresidue
end