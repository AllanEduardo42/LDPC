################################################################################
# Allan Eduardo Feitosa
# 01 Dec 2025
# VC-RBP

include("./RBP functions/findmaxedge.jl")
include("./RBP functions/calc_residue.jl")

function
    VC_RBP!(
        bitvector::Vector{Bool},
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},
        Lf::Vector{Float64},
        Nc::Vector{Vector{Int}},
        Nv::Vector{Vector{Int}},
        phi::Union{Vector{Float64},Nothing},
        num_reps::Int,
        alpha::Vector{Float64},
        Residues::Matrix{Float64},
        rbp_not_converged::Bool,
        msum_factor::Union{Float64,Nothing},
        # greediness::Vector{Int}
    )

    # greediness .= 0
    
    @fastmath for e in 1:num_reps

        # display("e = $e")

        # 1) Find largest residue  and coordenates
        cimax, vjmax = findmaxedge_VC(Residues,alpha,Nv)
        if vjmax == 0.0
            rbp_not_converged = false
            break # i.e., BP has converged
        end

        # greediness[vjmax] += 1

        # 2) Set to zero max residue
        Residues[cimax,vjmax] = 0.0

        Ncimax = Nc[cimax]        
        for vj in Ncimax
            if vj ≠ vjmax
                # 5) update messages Lr[cimax, vj]
                li = LinearIndices(Lr)[cimax,vj]
                alp = Residues[li]
                Lr[li] = calc_Lr(Ncimax,cimax,vj,Lq,msum_factor)
                # 6) calculate Residues
                Nvj = Nv[vj]
                Ld = calc_Ld(vj,Nvj,Lf,Lr)
                bitvector[vj] = signbit(Ld)
                for ci in Nvj
                    if ci ≠ cimax
                        li = LinearIndices(Lq)[ci,vj]
                        newlq = Ld - Lr[li]
                        oldlq = 2*atanh(Lq[li])
                        alp, residue = calc_residue_lq!(newlq,oldlq,alp)
                        Residues[li] = residue
                        Lq[li] = tanhLq(newlq,msum_factor)
                    end
                end
                alpha[vj] = alp
            end
        end
    end

    return rbp_not_converged
end

