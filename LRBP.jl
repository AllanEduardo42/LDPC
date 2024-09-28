
################################################################################
# Allan Eduardo Feitosa
# 26 set 2024
# local RBP Sum-Product Algorithm using min-sum to calculate the residues

include("minsum.jl")

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
        num_edges::Integer,
        Ldn::Vector{<:AbstractFloat},
        syndrome::Vector{Bool}
    )

    e = 1
    while e <= num_edges

        (cnmax,vnmax) = maxcoords
        Factors[cnmax,vnmax] *= pfactor
        Edges[cnmax,vnmax] += 1
        
        ### update Lr[cnmax,vnmax]
        update_Lr!(Lr,Lq,cnmax,vnmax,cn2vn)

        maxresidue = 0.0

        # update of Lq[vnmax,m], ∀m ≠ cnmax
        Ldn[vnmax] = Lf[vnmax]
        for m in vn2cn[vnmax]
            Ldn[vnmax] += Lr[m,vnmax]
            d[vnmax] = signbit(Ldn[vnmax])
        end
        for m in vn2cn[vnmax]
            if m ≠ cnmax
                Lq[vnmax,m] = Ldn[vnmax] - Lr[m,vnmax]
                maxresidue = minsum_LRBP!(
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

end
