################################################################################
# Allan Eduardo Feitosa
# 11 out 2024
# RBP Sum-Product Algorithm using min-sum to calculate the residues

include("calc_residues.jl")
include("findmaxcoords.jl")
include("update_list.jl")

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
        Lrn::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},        
        phi::Union{Vector{<:AbstractFloat},Nothing},
        Factors::Matrix{<:AbstractFloat},
        rbpfactor::AbstractFloat,
        num_edges::Integer,
        Ldn::Vector{<:AbstractFloat},
        samples::Union{Vector{<:Integer},Nothing},
        rng_sample::Union{AbstractRNG,Nothing},
        listsize1::Integer,
        listsize2::Integer,
        listres1::Union{Vector{<:AbstractFloat},Nothing},
        listm1::Union{Vector{<:Integer},Nothing},
        listn1::Union{Vector{<:Integer},Nothing},
        listres2::Union{Vector{<:AbstractFloat},Nothing},
        listm2::Union{Vector{<:Integer},Nothing},
        listn2::Union{Vector{<:Integer},Nothing},
        inlist::Union{Matrix{<:Integer},Nothing}
    )

    for e in 1:num_edges

        (cnmax,vnmax) = maxcoords

        if cnmax == 0
            cnmax = rand(rng_sample,1:length(cn2vn))
            vnmax = rand(rng_sample,cn2vn[cnmax])
        end

        _RBP_update_Lr!(cnmax,vnmax,Factors,rbpfactor,cn2vn,Lq,Lr,Ms,Lrn,signs)

        __RBP(Residues,cnmax,vnmax,listsize1,listres1,listm1,listn1,inlist)

        # update Ldn[vmax] and bitvector[vnmax]
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
                maxresidue = calc_residues!(alpha,Residues,maxcoords,0.0,
                                Factors,Ms,Lr,Lq,Lrn,signs,phi,vnmax,m,cn2vn,
                                listres1,listm1,listn1,listres2,listm2,listn2,
                                listsize1,listsize2,inlist)

            elseif listres1 !== nothing
                maxresidue = listres1[1]
                maxcoords[1], maxcoords[2] = listm1[1], listn1[1]
            end
        end

        @inbounds if listsize2 ≠ 0
            for k in 1:listsize2
                update_list!(inlist,listres1,listm1,listn1,listres2[k],
                    listm2[k],listn2[k],listsize1)
            end
            maxcoords[1], maxcoords[2] = listm1[1], listn1[1]
            listres2 .*= 0.0
            listm1 .*= 0
            listn1 .*= 0
        end

        maxresidue = findmaxcoords!(maxresidue,maxcoords,Residues,cn2vn,
            samples,rng_sample)

    end
end

# RBP and Random-RBP
function 
    __RBP(
        Residues::Matrix{<:AbstractFloat},
        cnmax::Integer,
        vnmax::Integer,
        ::Integer,
        ::Nothing,
        ::Nothing,
        ::Nothing,
        ::Nothing
    )

    @inbounds Residues[cnmax,vnmax] = 0.0
end

# Local-RBP
function 
    __RBP(
        ::Nothing,
        ::Integer,
        ::Integer,
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
        cnmax::Integer,
        vnmax::Integer,
        listsize::Integer,
        listres1::Vector{<:AbstractFloat},
        listm1::Vector{<:Integer},
        listn1::Vector{<:Integer},
        inlist::Matrix{Bool}
    )

    @inbounds inlist[cnmax,vnmax] = false
    @inbounds for i in 1:listsize
        listres1[i] = listres1[i+1]
        listm1[i] = listm1[i+1]
        listn1[i] = listn1[i+1]
    end    
end

function 
    _RBP_update_Lr!(
        cnmax::Integer,
        vnmax::Integer,
        Factors::Matrix{<:AbstractFloat},
        rbpfactor::AbstractFloat,
        ::Vector{Vector{T}} where {T<:Integer},
        ::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Ms::Matrix{<:AbstractFloat},
        ::Vector{<:AbstractFloat},
        ::Nothing
    )

    @fastmath @inbounds Factors[cnmax,vnmax] *= rbpfactor

    @inbounds Lr[cnmax,vnmax] = Ms[cnmax,vnmax]

end

function 
    _RBP_update_Lr!(
        cnmax::Integer,
        vnmax::Integer,
        Factors::Matrix{<:AbstractFloat},
        rbpfactor::AbstractFloat,
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        ::Matrix{<:AbstractFloat},
        ::Union{Vector{<:AbstractFloat},Nothing},
        ::Vector{Bool}
    )

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