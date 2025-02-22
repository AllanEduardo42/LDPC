################################################################################
# Allan Eduardo Feitosa
# 22 Feb 2024
# RBP Sum-Product Algorithm with residual decay factor

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/calc_residues.jl")
include("./RBP functions/findmaxedge.jl")

function
    RBP!(
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
        Ldn::Vector{<:AbstractFloat},
        Ms::Matrix{<:AbstractFloat},
        Factors::Matrix{<:AbstractFloat},        
        all_max_res_alt::Union{Vector{<:AbstractFloat},Nothing},
        test::Bool,
        address::Union{Matrix{<:Integer},Nothing},
        addressinv::Union{Matrix{<:Integer},Nothing},
        residues::Union{Vector{<:AbstractFloat},Nothing},
    )

    @inbounds @fastmath for e in 1:num_edges

        # display("e = $e")

        # 1) Find maximum residue
        max_edge = findmaxedge(residues)
        if max_edge == 0
            break
        else
            cnmax = address[1,max_edge]
            vnmax = address[2,max_edge]
        end

        # 2) for optimization, get the linear indice of the maximum residue
        lmax = LinearIndices(Factors)[cnmax,vnmax]

        # 2.5) if in test mode, store the values of the maximum residues
        if test
            all_max_res_alt[e] = residues[addressinv[lmax]]   
        end 

        # 3) Decay the RBP factor corresponding to the maximum residue
        Factors[lmax] *= decayfactor

        # 4) update check to node message Lr[cnmax,vnmax]
        RBP_update_Lr!(lmax,cnmax,vnmax,cn2vn,Lq,Lr,Ms,Lrn,signs,phi)

        # 5) set maximum residue to zero or remove it from the list
        residues[addressinv[lmax]] = 0.0

        # 6) update Ldn[vmax] and bitvector[vnmax]
        Ldn[vnmax], nl = calc_Ld(vnmax,vn2cn,Lf[vnmax],Lr)
        bitvector[vnmax] = signbit(Ldn[vnmax])

        # 7) update vn2cn messages Lq[vnmax,m], ∀m ≠ cnmax, and calculate residues
        leaf = true # suppose node vnmax is a leaf in the graph
        for m in vn2cn[vnmax]
            if m ≠ cnmax
                leaf = false # vnmax is not a leaf
                Lq[vnmax,m] = Ldn[vnmax] - Lr[nl+m]
                calc_residues!(addressinv,residues,Factors,Ms,Lr,Lq,Lrn,signs,
                               phi,vnmax,m,cn2vn)
            end
        end

        # 8) if vnmax is a leaf in the graph
        if leaf
            Lq[vnmax,cnmax] = Ldn[vnmax] - Lr[nl+cnmax]
            calc_residues!(addressinv,residues,Factors,Ms,Lr,Lq,Lrn,signs,
                            phi,vnmax,m,cn2vn)
        end
    end
end

