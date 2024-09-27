
################################################################################
# Allan Eduardo Feitosa
# 26 set 2024
# local RBP Sum-Product Algorithm using min-sum to calculate the residues

include("min_sum.jl")

function
    LRBP!(
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
        num_edges::Integer
    )

    e = 1
    while e <= num_edges

        (cnmax,vnmax) = maxcoords
        Factors[cnmax,vnmax] *= pfactor
        Edges[cnmax,vnmax] += 1
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
        maxresidue = 0.0
        
        checknodes_vnmax = vn2cn[vnmax]
        for m in checknodes_vnmax
            if m ≠ cnmax
                # vertical update of Lq[vnmax,m], ∀cn ≠ cnmax
                @inbounds Lq[vnmax,m] = Lf[vnmax]
                for m2 in checknodes_vnmax
                    if m2 ≠ m
                        @inbounds @fastmath Lq[vnmax,m] += Lr[m2,vnmax]
                    end
                end
                maxresidue = min_sum_lRBP!(
                    maxcoords,
                    maxresidue,
                    Factors,
                    Lr,
                    Lq,
                    signs,
                    vnmax,
                    m,
                    cn2vn
                )
            end
        end

        if maxresidue == 0.0
            break
        end

        e += 1

    end

    MAP!(d,vn2cn,Lf,Lr)

end
