################################################################################
# Allan Eduardo Feitosa
# 31 Mar 2025
# RBP Sum-Product Algorithm with residual decaying factor

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/findmaxnode.jl")

function
    VN_RBP!(
        bitvector::Vector{Bool},
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        Lrn::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{<:AbstractFloat},Nothing},
        decayfactor::AbstractFloat,
        num_edges::Integer,
        Ms::Matrix{<:AbstractFloat},
        Factors::Vector{<:AbstractFloat},
        residues::Vector{<:AbstractFloat},
        relative::Bool,
        rbp_not_converged::Bool
    )

    @fastmath @inbounds for e = 1:num_edges

        # display("e = $e")

        # display(sort(residues,rev=true))

        vnmax = findmaxnode(residues)
        if vnmax == 0
            rbp_not_converged  = false
            break # i.e., RBP has converged
        end
        residues[vnmax] = 0.0
        Factors[vnmax] *= decayfactor

        ln = LinearIndices(Ms)[1,vnmax]-1
        for m in vn2cn[vnmax]
            Lr[ln+m] = Ms[ln+m]
        end        

        bitvector[vnmax] = update_Lq!(Lq,Lr,Lf[vnmax],vnmax,vn2cn[vnmax],Lrn)

        for m in vn2cn[vnmax] 
            # calculate the new check to node messages
            vns = cn2vn[m]
            update_Lr!(Ms,Lq,m,vns,Lrn,signs,phi)
            for n in vns
                if n ≠ vnmax
                    li = LinearIndices(Ms)[m,n]
                    residue = Ms[li] - Lr[li]
                    if relative   
                        old = Lr[li] + Lq[li]  
                        residue /= old
                    end
                    residues[n] = abs(residue)*Factors[n]
                end
            end
        end
    end

    return rbp_not_converged

end

