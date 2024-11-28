################################################################################
# Allan Eduardo Feitosa
# 11 out 2024
# RBP Sum-Product Algorithm using min-sum to calculate the residues

include("RBP_functions.jl")
include("findmaxcoords.jl")

#RBP
function
    RBP!(
        alpha::AbstractFloat,
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
        listsize1::Integer,
        listsize2::Integer,
        listres1::Union{Vector{<:AbstractFloat},Nothing},
        listadd1::Union{Matrix{<:Integer},Nothing},
        listres2::Union{Vector{<:AbstractFloat},Nothing},
        listadd2::Union{Matrix{<:Integer},Nothing},
        listaddinv1::Union{Matrix{<:Integer},Nothing},
        inlist1::Union{Matrix{<:Integer},Nothing}
    )

    for e in 1:num_edges

        # display(listadd1)
        # display(listadd2)

        if iszero(maxcoords) || maxresidue == 0.0
            maxcoords[1] = rand(rng_sample,1:length(cn2vn))
            maxcoords[2] = rand(rng_sample,cn2vn[maxcoords[1]])
        end

        _RBP_update_Lr!(maxcoords,Factors,rbpfactor,cn2vn,Lq,Lr)

        __RBP(Residues,maxcoords,listsize1,listres1,listadd1,listaddinv1,inlist1)

        maxresidue = _RBP_update_vn2cn!(alpha,Residues,maxcoords,0.0,Factors,Lf,Ldn,
            bitvector,vn2cn,cn2vn,Ms,Lr,Lq,signs,listsize1,listsize2,listres1,
            listadd1,listres2,listadd2,listaddinv1,inlist1)

        maxresidue = findmaxcoords!(maxresidue,maxcoords,Residues,cn2vn,
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
        listres1::Vector{<:AbstractFloat},
        listadd1::Matrix{<:Integer},
        listaddinv::Matrix{<:Integer},
        inlist1::Matrix{Bool}
    )

    # println()
    # println("maxcoords:")
    # display(maxcoords)
    # println()

    @inbounds inlist1[maxcoords[1],maxcoords[2]] = false
    @inbounds listaddinv[maxcoords[1],maxcoords[2]] = 0
    @inbounds for i in 1:listsize
        listres1[i] = listres1[i+1]
        m = listadd1[1,i+1]
        n = listadd1[2,i+1]
        listadd1[1,i] = m
        listadd1[2,i] = n
        if listadd1[1,i] != 0
            listaddinv[m,n] = i
        end
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
    @fastmath @inbounds for n in cn2vn[cnmax]
        if n != vnmax
            pLr *= tanh(0.5*Lq[n,cnmax])
        end
    end    
    if @fastmath abs(pLr) < 1 
        @fastmath @inbounds Lr[cnmax,vnmax] = 2*atanh(pLr)
    elseif pLr > 0
        @fastmath @inbounds Lr[cnmax,vnmax] = INFFLOAT
    else
        @fastmath @inbounds Lr[cnmax,vnmax] = NINFFLOAT
    end

end

function 
    _RBP_update_vn2cn!(
        alpha::AbstractFloat,
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
        listsize1::Integer,
        listsize2::Integer,
        listres1::Union{Vector{<:AbstractFloat},Nothing},
        listadd1::Union{Matrix{<:Integer},Nothing},
        listres2::Union{Vector{<:AbstractFloat},Nothing},
        listadd2::Union{Matrix{<:Integer},Nothing},
        listaddinv1::Union{Matrix{<:Integer},Nothing},
        inlist1::Union{Matrix{<:Integer},Nothing}       
    )

    # update Ldn[vmax] and bitvector[vnmax]
    (cnmax,vnmax) = maxcoords
    @inbounds Ldn[vnmax] = Lf[vnmax]
    @fastmath @inbounds for m in vn2cn[vnmax]
        Ldn[vnmax] += Lr[m,vnmax]
        bitvector[vnmax] = signbit(Ldn[vnmax])
    end

    @fastmath @inbounds for m in vn2cn[vnmax]
        # if m ≠ cnmax || length(vn2cn[vnmax]) == 1
        if m ≠ cnmax
            # update vn2cn messages Lq[vnmax,m], ∀m ≠ cnmax
            Lq[vnmax,m] = Ldn[vnmax] - Lr[m,vnmax]
            maxresidue = calc_residues!(alpha,Residues,maxcoords,maxresidue,Factors,
                                        Ms,Lr,Lq,signs,vnmax,m,cn2vn,listres1,
                                        listadd1,listres2,listadd2,listaddinv1,
                                        listsize1,listsize2,inlist1)
        elseif listres1 !== nothing
            maxresidue = listres1[1]
            maxcoords[1] = listadd1[1,1]
            maxcoords[2] = listadd1[2,1]
        end
    end

    # println()
    # println("### lists ###")
    # println()
    # println("list 1")
    # display(listres1)
    # display(listadd1)
    # println()
    # println("list 2")
    # display(listres2)
    # display(listadd2)

    if listsize2 != 0
        @fastmath @inbounds for k=1:listsize2
            y = listres2[k]
            # for i=1:listsize1
            if y > listres1[listsize1]
                if y ≥ listres1[1]    
                    i = 1
                else
                    d = listsize1 >> 1
                    i = d
                    while d > 1
                        d >>= 1
                        if y ≥ listres1[i]
                            i -= d
                        else
                            i += d
                        end
                    end
                    if y < listres1[i]
                        i += 1
                    end
                end
                # println("i = $i")
                mm = listadd1[1,end-1]
                nn = listadd1[2,end-1]
                if mm ≠ 0
                    inlist1[mm,nn] = false
                    listaddinv1[mm,nn] = 0
                end
                for j=listsize1:-1:i+1
                    listres1[j] = listres1[j-1]
                    listadd1[1,j] = listadd1[1,j-1]
                    listadd1[2,j] = listadd1[2,j-1]
                    if listadd1[1,j] != 0
                        listaddinv1[listadd1[1,j],listadd1[2,j]] = j
                    end
                end
                listadd1[1,i] = listadd2[1,k]
                listadd1[2,i] = listadd2[2,k]
                listaddinv1[listadd2[1,k],listadd2[2,k]] = i
                listres1[i] = y
                inlist1[listadd2[1,k],listadd2[2,k]] = true
            end
        end
        @inbounds maxcoords[1] = listadd1[1,1]
        @inbounds maxcoords[2] = listadd1[2,1]
        listres2 .*= 0.0
        listadd2 .*= 0
    end

    # println()
    # println("list 1 updated")
    # display(listres1)
    # display(listadd1)

    return maxresidue
end