################################################################################
# Allan Eduardo Feitosa
# 31 Mar 2025
# RBP Sum-Product Algorithm with residual decaying factor

include("./RBP functions/findmaxnode.jl")
include("./RBP functions/calc_residue.jl")

function
    VN_RBP!(
        bitvector::Vector{Bool},
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        Nc::Vector{Vector{T}} where {T<:Integer},
        Nv::Vector{Vector{T}} where {T<:Integer},
        aux::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{<:AbstractFloat},Nothing},
        decayfactor::AbstractFloat,
        num_edges::Integer,
        newLr::Matrix{<:AbstractFloat},
        Factors::Vector{<:AbstractFloat},
        residues::Vector{<:AbstractFloat},
        bp_not_converged::Bool
    )

    @fastmath @inbounds for e in 1:num_edges

        # display("e = $e")

        vjmax = findmaxnode(residues)
        if vjmax == 0
            bp_not_converged = false
            break # i.e., BP has converged
        end
        Factors[vjmax] *= decayfactor
        residues[vjmax] = 0.0

        for ci in Nv[vjmax]
            li = LinearIndices(Lr)[ci,vjmax]
            RBP_update_Lr!(li,Lr,newLr,ci,vjmax,Nc[ci],Lq,aux,signs,phi)
        end   
        
        Nvjmax = Nv[vjmax]
        Ld = calc_Ld(vjmax,Nvjmax,Lf,Lr)
        bitvector[vjmax] = signbit(Ld)

        for ci in Nv[vjmax]
            # update Nv messages Lq[ci,vjmax]
            li = LinearIndices(Lq)[ci,vjmax]
            Lq[li] = Ld - Lr[li]
            # calculate the new check to node messages
            Nci = Nc[ci]
            A, B, C, D = calc_ABCD!(aux,signs,phi,Lq,ci,Nci)
            for vj in Nci
                if vj ≠ vjmax
                    newlr = calc_Lr(A,B,C,D,vj,aux,signs,phi)
                    li = LinearIndices(newLr)[ci,vj] 
                    newLr[li] = newlr   
                    residues[vj] = calc_residue(newlr,Lr[li],Factors[vj])
                end
            end
        end
    end

    return bp_not_converged

end

function
    VN_RBP!(
        bitvector::Vector{Bool},
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        Nc::Vector{Vector{T}} where {T<:Integer},
        Nv::Vector{Vector{T}} where {T<:Integer},
        ::Nothing,
        ::Nothing,
        ::Nothing,
        decayfactor::AbstractFloat,
        num_edges::Integer,
        newLr::Matrix{<:AbstractFloat},
        Factors::Vector{<:AbstractFloat},
        residues::Vector{<:AbstractFloat},
        bp_not_converged::Bool
    )

    @fastmath @inbounds for e in 1:num_edges

        # display("e = $e")

        vjmax = findmaxnode(residues)
        if vjmax == 0
            bp_not_converged = false
            break # i.e., BP has converged
        end
        Factors[vjmax] *= decayfactor
        residues[vjmax] = 0.0

        for ci in Nv[vjmax]
            li = LinearIndices(Lr)[ci,vjmax]
            Lr[li] = newLr[li]
        end 

        Nvjmax = Nv[vjmax]
        for ci in Nvjmax
            # 5) update Nv messages Lq[ci,vnmax]
            Lq[ci,vjmax] = calc_Lq(Nvjmax,ci,vjmax,Lr,Lf)
            # 6) calculate residues
            Nci = Nc[ci]
            for vj in Nci
                if vj ≠ vjmax
                    newlr = calc_Lr(Nci,ci,vj,Lq)
                    li = LinearIndices(Lr)[ci,vj]
                    newLr[li] = newlr
                    residues[vj] = calc_residue_raw(newlr,Lr[li],Factors[vj])
                end
            end
        end
    end

    # 7) update bitvector
    for vj in eachindex(Nv)
        ci = Nv[vj][1]
        li = LinearIndices(Lr)[ci,vj]
        bitvector[vj] = signbit(Lr[li] + Lq[li])
    end

    return bp_not_converged
end

