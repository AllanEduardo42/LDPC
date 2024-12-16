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
        address::Union{Matrix{<:Integer},Nothing},
        addressinv::Union{Matrix{<:Integer},Nothing},
        residues::Union{Vector{<:AbstractFloat},Nothing},
        bitvector::Vector{Bool},
        Lr::Matrix{<:AbstractFloat},
        Ms::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        Lrn::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},        
        phi::Union{Vector{<:AbstractFloat},Nothing},
        Factors::Matrix{<:AbstractFloat},
        decayfactor::AbstractFloat,
        num_edges::Integer,
        Ldn::Vector{<:AbstractFloat},
        rng_sample::AbstractRNG,
        listsize1::Integer,
        listsize2::Integer,
        listres1::Vector{<:AbstractFloat},
        listm1::Vector{<:Integer},
        listn1::Vector{<:Integer},
        listres2::Union{Vector{<:AbstractFloat},Nothing},
        listm2::Union{Vector{<:Integer},Nothing},
        listn2::Union{Vector{<:Integer},Nothing},
        inlist::Union{Matrix{<:Integer},Nothing}
    )

    @fastmath @inbounds for e in 1:num_edges

        # println(listres1[1])

        if listres1[1] == 0
            break
        end
        
        cnmax = listm1[1]
        vnmax = listn1[1]

        if cnmax == 0 || listm1[1]== -1
            cnmax = rand(rng_sample,1:length(cn2vn))
            vnmax = rand(rng_sample,cn2vn[cnmax])
        end

        lmax = LinearIndices(Factors)[cnmax,vnmax]

        Factors[lmax] *= decayfactor

        ### update Lr[cnmax,vnmax]
        _RBP_update_Lr!(lmax,cnmax,vnmax,cn2vn,Lq,Lr,Ms,Lrn,signs)

        __RBP(addressinv,residues,lmax,listsize1,listres1,listm1,listn1,inlist)

        # update Ldn[vmax] and bitvector[vnmax]
        Ldn[vnmax] = Lf[vnmax]
        for m in vn2cn[vnmax]
            Ldn[vnmax] += Lr[m,vnmax]
            bitvector[vnmax] = signbit(Ldn[vnmax])
        end

        leaf = true
        for m in vn2cn[vnmax]
            if m ≠ cnmax
                leaf = false
                # update vn2cn messages Lq[vnmax,m], ∀m ≠ cnmax
                Lq[vnmax,m] = Ldn[vnmax] - Lr[m,vnmax]
                calc_residues!(addressinv,residues,Factors,Ms,Lr,Lq,Lrn,signs,phi,vnmax,m,
                    cn2vn,listres1,listm1,listn1,listres2,listm2,listn2,
                    listsize1,listsize2,inlist)
            end
        end

        if leaf
            if listsize1 == 1
                listm1[1] = 0
            end
        end

        if listsize2 ≠ 0
            for k in listsize2:-1:1
            # for k in 1:listsize2
                update_list!(inlist,listres1,listm1,listn1,listres2[k],
                    listm2[k],listn2[k],listsize1)
            end
            listres2 .*= 0.0
            listm2 .*= 0
            listn2 .*= 0
        end

        findmaxcoords!(address,residues,listres1,listm1,listn1)

    end
end

# RBP
function 
    __RBP(
        addressinv::Matrix{<:Integer},
        residues::Vector{<:AbstractFloat},
        lmax::Integer,
        ::Integer,
        listres1::Vector{<:AbstractFloat},
        ::Vector{<:Integer},
        ::Vector{<:Integer},
        ::Nothing
    )

    @inbounds listres1[1] = 0
    @inbounds residues[addressinv[lmax]] = 0.0

end

# Local-RBP
function 
    __RBP(
        ::Nothing,
        ::Nothing,
        ::Integer,
        ::Integer,
        listres1::Vector{<:AbstractFloat},
        ::Vector{<:Integer},
        ::Vector{<:Integer},
        ::Nothing
    )
    
    @inbounds listres1[1] = 0

end

# List-RBP
function 
    __RBP(
        ::Nothing,
        ::Nothing,
        lmax::Integer,
        listsize::Integer,
        listres1::Vector{<:AbstractFloat},
        listm1::Vector{<:Integer},
        listn1::Vector{<:Integer},
        inlist::Matrix{Bool}
    )

    @inbounds inlist[lmax] = false
    @inbounds for i in 1:listsize
        listres1[i] = listres1[i+1]
        listm1[i] = listm1[i+1]
        listn1[i] = listn1[i+1]
    end    
end

function 
    _RBP_update_Lr!(
        lmax::Integer,
        ::Integer,
        ::Integer,
        ::Vector{Vector{T}} where {T<:Integer},
        ::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Ms::Matrix{<:AbstractFloat},
        ::Vector{<:AbstractFloat},
        ::Nothing
    )

    @inbounds Lr[lmax] = Ms[lmax]

end

function 
    _RBP_update_Lr!(
        lmax::Integer,
        cnmax::Integer,
        vnmax::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        ::Matrix{<:AbstractFloat},
        ::Union{Vector{<:AbstractFloat},Nothing},
        ::Vector{Bool}
    )

    pLr = 1.0
    @fastmath @inbounds for n in cn2vn[cnmax]
        if n != vnmax
            pLr *= tanh(0.5*Lq[n,cnmax])
        end
    end    
    if @fastmath abs(pLr) < 1 
        @fastmath @inbounds Lr[lmax] = 2*atanh(pLr)
    elseif pLr > 0
        @fastmath @inbounds Lr[lmax] = INFFLOAT
    else
        @fastmath @inbounds Lr[lmax] = NINFFLOAT
    end

end