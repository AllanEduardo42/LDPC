
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
        Edges::Matrix{<:Integer},
        Factors::Matrix{<:AbstractFloat},
        pfactor::AbstractFloat,
        num_edges::Integer,
        Residues::Matrix{<:AbstractFloat},
        samples::Vector{<:Integer}
    )

    e = 1
    while e <= num_edges

        maxresidue = find_maxresidue_coords!(
            maxcoords,
            Residues,
            cn2vn,
            samples)

        if maxresidue == 0.0 # if RBP has converged
            break
        end

        (cnmax,vnmax) = maxcoords
        Factors[cnmax,vnmax] *= pfactor
        Edges[cnmax,vnmax] += 1
        ### update Lr[cnmax,vnmax]
        update_Lr!(Lr,Lq,cnmax,vnmax,cn2vn)

        # we don't use the m-to-n message correspoding to (cnmax,vnmax) anymore.
        # Thus we make the largest residue equal to zero:
        Residues[cnmax,vnmax] = 0.0
        
        checknodes_vnmax = vn2cn[vnmax]
        for m in checknodes_vnmax
            if m ≠ cnmax
                # update Lq[vnmax,m], ∀cn ≠ cnmax
                update_Lq!(Lq,Lr,Lf,m,vnmax,vn2cn)
                # if any new residue estimate is larger than the previously estimated maximum 
                # residue than update the value of maxresidue and maxcoords.
                minsum_RBP!(Residues,Factors,Lr,Lq,signs,vnmax,m,cn2vn)
            end
        end

        e +=1

    end

    MAP!(d,vn2cn,Lf,Lr)

end
function 
    update_Lr!(
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        cnmax::Integer,
        vnmax::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer}
    )

    pLr = 1.0
    for n in cn2vn[cnmax]
        if n != vnmax
            @inbounds @fastmath pLr *= tanh(0.5*Lq[n,cnmax])
        end
    end    
    if abs(pLr) < 1 
        @inbounds @fastmath Lr[cnmax,vnmax] = 2*atanh(pLr)
    end
end

function 
    update_Lq!(
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        m::Integer,
        vnmax::Integer,
        vn2cn::Vector{Vector{T}} where {T<:Integer}
    )

    @inbounds Lq[vnmax,m] = Lf[vnmax]
    for m2 in vn2cn[vnmax]
        if m2 ≠ m
            @inbounds @fastmath Lq[vnmax,m] += Lr[m2,vnmax]
        end
    end
end

function
    find_maxresidue_coords!(
        maxcoords::Vector{<:Integer},
        Residues::Matrix{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        samples::Vector{<:Integer}
    )

    maxresidue = 0
    rand!(samples,1:size(Residues,1))
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