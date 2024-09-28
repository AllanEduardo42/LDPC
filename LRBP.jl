
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
        num_edges::Integer
    )

    e = 1
    while e <= num_edges

        (cnmax,vnmax) = maxcoords
        Factors[cnmax,vnmax] *= pfactor
        Edges[cnmax,vnmax] += 1
        ### update Lr[cnmax,vnmax]
        update_Lr!(Lr,Lq,cnmax,vnmax,cn2vn)

        maxresidue = 0.0
        
        for m in vn2cn[vnmax]
            if m ≠ cnmax
                # update of Lq[vnmax,m], ∀cn ≠ cnmax
                update_Lq!(Lq,Lr,Lf,m,vnmax,vn2cn)
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

    MAP!(d,vn2cn,Lf,Lr)

end
