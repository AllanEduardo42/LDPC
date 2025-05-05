################################################################################
# Allan Eduardo Feitosa
# 07 Apr 2025
# RBP and List-RBP Algorithm with residual decaying factor

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/findmaxedge.jl")
include("./RBP functions/update_lists.jl")
include("./RBP functions/remove_residue!.jl")

function
    RBP_raw!(
        bitvector::Vector{Bool},
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        ::Union{Vector{<:AbstractFloat},Nothing},
        ::Union{Vector{Bool},Nothing},
        ::Union{Vector{<:AbstractFloat},Nothing},
        decayfactor::AbstractFloat,
        num_edges::Integer,
        Ms::Matrix{<:AbstractFloat},
        Factors::Matrix{<:AbstractFloat},
        coords::Matrix{<:Integer},
        indices::Matrix{<:Integer},
        residues::Vector{<:AbstractFloat},
        local_residues::Union{Vector{<:AbstractFloat},Nothing},
        ::Union{Matrix{<:Integer},Nothing},
        ::Vector{<:Integer},
        ::Bool,
        bp_not_converged::Bool
    )

    @fastmath @inbounds for e in 1:num_edges

        # display("e = $e")

        # 1) Find largest residue  and coordenates
        max_edge, max_residue = findmaxedge(residues,local_residues)
        if max_residue == 0.0
            if max_edge == 0
                bp_not_converged = false
                break # i.e., RBP has converged
            end
        else
            cnmax = coords[1,max_edge]
            vnmax = coords[2,max_edge]
            limax = coords[3,max_edge]
        end

        # 2) Decay the RBP factor corresponding to the maximum residue
        Factors[limax] *= decayfactor
        # 3) update check to node message Lr[cnmax,vnmax]
        pLr = 1.0
        for n2 in cn2vn[cnmax]  
            if n2 ≠ vnmax
                pLr *= tanh(0.5*Lq[cnmax,n2])
            end
        end
        Lr[cnmax,vnmax] = 2*atanh(pLr)
        # 4) set maximum residue to zero
        residues[max_edge] = 0.0

        # 5) update vn2cn messages Lq[vnmax,m] and bitvector[vnmax]
        cns = vn2cn[vnmax]
        Ld = calc_Ld(vnmax,cns,Lf[vnmax],Lr)
        bitvector[vnmax] = signbit(Ld)

        # 6) calculate residues
        for m in cns
            if m ≠ cnmax
                Lq[m,vnmax] = Ld - Lr[m,vnmax]
                vns = cn2vn[m]    
                # calculate the residues
                for n in vns
                    if n ≠ vnmax
                        pLr = 1.0
                        for n2 in vns
                            if n2 ≠ n
                                pLr *= tanh(0.5*Lq[m,n2])
                            end
                        end
                        residues[indices[m,n]] = calc_residue(pLr,Lr,Factors,m,n)
                    end
                end
            end
        end
    end

    return bp_not_converged
end

function calc_residue(
    pLr::AbstractFloat,
    Lr::Matrix{<:AbstractFloat},
    Factors::Matrix{<:AbstractFloat},
    m::Integer,
    n::Integer
)

    @inbounds begin
        x = 2*atanh(pLr) - Lr[m,n]
        if isnan(x)
            return 0.0
        else
            @fastmath return abs(x)*Factors[m,n]
        end
    end
end