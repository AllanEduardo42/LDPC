
################################################################################
# Allan Eduardo Feitosa
# 23 set 2024
# RBP Sum-Product Algorithm using min-sum to calculate the residues

include("minsum_RBP.jl")

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
        factor::AbstractFloat,
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
        Factors[cnmax,vnmax] *= factor

        ### update Lr[cnmax,vnmax]
        pLr = 1.0
        for n in cn2vn[cnmax]
            if n != vnmax
                @inbounds @fastmath pLr *= tanh(0.5*Lq[n,cnmax])
            end
        end    
        if abs(pLr) < 1 
            @inbounds @fastmath Lr[cnmax,vnmax] = 2*atanh(pLr)
        end

        # we don't use the m-to-n message correspoding to (cnmax,vnmax) anymore.
        # Thus we make the largest residue equal to zero:
        if lrbp
            maxresidue = 0.0
        else
            Residues[cnmax,vnmax] = 0.0
        end
        
        # update Ldn[vmax] and d[vnmax]
        Ldn[vnmax] = Lf[vnmax]
        for m in vn2cn[vnmax]
            Ldn[vnmax] += Lr[m,vnmax]
            d[vnmax] = signbit(Ldn[vnmax])
        end

        for m in vn2cn[vnmax]
            if m ≠ cnmax
                # update vn2cn messages Lq[vnmax,m], ∀m ≠ cnmax
                Lq[vnmax,m] = Ldn[vnmax] - Lr[m,vnmax]
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

function
    find_maxresidue_coords!(
        maxcoords::Vector{<:Integer},
        Residues::Matrix{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        samples::Vector{<:Integer},
        rng_sample::AbstractRNG
    )

    maxresidue = 0
    M = size(Residues,1)
    rand!(rng_sample,samples,1:M)
    for m in samples
        for n in cn2vn[m]
            if Residues[m,n] > maxresidue
                maxresidue = Residues[m,n]
                maxcoords[1] = m
                maxcoords[2] = n
            end
        end
    end

    return maxresidue
end

function
    find_maxresidue_coords!(
        maxcoords::Vector{<:Integer},
        Residues::Matrix{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        samples::Nothing,
        rng_sample::AbstractRNG
    )

    maxresidue = 0
    for m in eachindex(cn2vn)
        for n in cn2vn[m]
            if Residues[m,n] > maxresidue
                maxresidue = Residues[m,n]
                maxcoords[1] = m
                maxcoords[2] = n
            end
        end
    end

    return maxresidue
end