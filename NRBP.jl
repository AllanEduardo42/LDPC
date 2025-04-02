################################################################################
# Allan Eduardo Feitosa
# 31 Mar 2025
# RBP Sum-Product Algorithm with residual decaying factor

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/calc_local_residues.jl")
include("./RBP functions/findmaxnode.jl")
include("./RBP functions/decay.jl")

function
    NRBP!(
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
        all_max_res_alt::Union{Vector{<:AbstractFloat},Nothing},
        test::Bool,
        residues::Vector{<:AbstractFloat}
    )

    @fastmath @inbounds for e = 1:num_edges

        # display("e = $e")

        # display(sort(residues,rev=true))

        vnmax, maxresidue = findmaxnode(residues)
        if vnmax == 0
            break # i.e., RBP has converged
        end
        residues[vnmax] = 0.0
        Factors[vnmax] *= decayfactor

        ln = LinearIndices(Ms)[1,vnmax]-1
        for m in vn2cn[vnmax]
            Lr[ln+m] = Ms[ln+m]
        end        

        bitvector[vnmax] = update_Lq!(Lq,Lr,Lf[vnmax],vnmax,vn2cn,Lrn)

        for m in vn2cn[vnmax] 
            # calculate the new check to node messages
            update_Lr!(Ms,Lq,m,cn2vn,Lrn,signs,phi)
            for n in cn2vn[m]
                if n ≠ vnmax
                    li = LinearIndices(Ms)[m,n]      
                    old = Lr[li] + Lq'[li]  
                    residues[n] = abs((Ms[li] - Lr[li])/old)*Factors[n]
                end
            end
        end
    end
end

