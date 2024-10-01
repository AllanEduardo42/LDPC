
################################################################################
# Allan Eduardo Feitosa
# 23 set 2024
# RBP Sum-Product Algorithm using min-sum to calculate the residues

include("RBP_functions.jl")

function
    RBP!(
        d::Vector{Bool},
        Lr::Matrix{<:AbstractFloat},
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
        Residues::Union{Matrix{<:AbstractFloat},Nothing},
        samples::Union{Vector{<:Integer},Nothing},
        rng_sample::AbstractRNG,
        lrbp::Bool
    )

    e = 1
    while e <= num_edges

        if !lrbp

            maxresidue = find_maxresidue_coords!(
                maxcoords,
                Residues,
                cn2vn,
                samples,
                rng_sample)

            if maxresidue == 0.0 # if RBP has converged
                break
            end
        end

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

        # we don't use the m-to-n message correspoding to (cnmax,vnmax) anymore.
        # Thus we make the largest residue equal to zero:
        if lrbp
            maxresidue = 0.0
        else
            @inbounds Residues[cnmax,vnmax] = 0.0
        end
        
        # update Ldn[vmax] and d[vnmax]
        @inbounds Ldn[vnmax] = Lf[vnmax]
        for m in vn2cn[vnmax]
            @fastmath @inbounds Ldn[vnmax] += Lr[m,vnmax]
            @fastmath @inbounds d[vnmax] = signbit(Ldn[vnmax])
        end

        for m in vn2cn[vnmax]
            if m ≠ cnmax
                # update vn2cn messages Lq[vnmax,m], ∀m ≠ cnmax
                @fastmath @inbounds Lq[vnmax,m] = Ldn[vnmax] - Lr[m,vnmax]
                # if any new residue estimate is larger than the previously estimated maximum 
                # residue than update the value of maxresidue and maxcoords.
                maxresidue = minsum_RBP!(
                Residues,
                maxcoords,
                maxresidue,
                Factors,
                Lr,
                Lq,
                signs,
                vnmax,
                m,
                cn2vn)
            end
        end

        if maxresidue == 0.0 #breaks the loop in LRBP mode
            break
        end

        e +=1

    end

end