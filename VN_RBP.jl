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
        Lrj::Vector{<:AbstractFloat},
        decayfactor::AbstractFloat,
        num_edges::Integer,
        newLr::Matrix{<:AbstractFloat},
        Factors::Vector{<:AbstractFloat},
        residues::Vector{<:AbstractFloat},
        bp_not_converged::Bool
    )

    # count = 0

    @fastmath @inbounds for e in 1:num_edges
    # @fastmath @inbounds while count < 200000

        # display("e = $e")

        vjmax = findmaxnode(residues)
        if vjmax == 0
            bp_not_converged = false
            break # i.e., BP has converged
        end
        Factors[vjmax] *= decayfactor
        residues[vjmax] = 0.0

        for ci in Nv[vjmax]
            li = LinearIndices(newLr)[ci,vjmax]
            Lr[li] = newLr[li]
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
            pLr, countzeros, vj0 = calc_pLr(Lq,ci,Nci,Lrj) 
            for vj in Nci
                if vj ≠ vjmax
                    # count += 1
                    li = LinearIndices(newLr)[ci,vj]
                    newlr = fast_Lr(Lrj,pLr,countzeros,vj0,vj)
                    newLr[li]= newlr   
                    residue = newlr- Lr[li]
                    residues[vj] = abs(residue)*Factors[vj]
                end
            end
        end
    end

    # display(count)

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
                    li = LinearIndices(Lr)[ci,vj]
                    oldLr = Lr[li]
                    newlr = calc_Lr(Nci,ci,vj,Lq)
                    newLr[li] = newlr
                    residues[vj] = calc_residue(newlr,oldLr,Factors[vj])
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

