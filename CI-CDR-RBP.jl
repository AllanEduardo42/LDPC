################################################################################
# Allan Eduardo Feitosa
# 04 Mar 2026
# CI-CDR-RBP Algorithm with residual decaying factor

include("./RBP functions/findmaxedge.jl")
include("./RBP functions/calc_residue.jl")

function
    CI_CDR_RBP!(
        bitvector::Vector{Bool},
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},
        Lf::Vector{Float64},
        Nc::Vector{Vector{Int}},
        Nv::Vector{Vector{Int}},
        phi::Union{Vector{Float64},Nothing},
        msum_factor::Union{Float64,Nothing},
        msum2::Bool,
        num_reps::Int,
        newLr::Matrix{Float64},
        Residues::Matrix{Float64},
        alpha::Vector{Float64},
        rbp_not_converged::Bool,
        Dn::Vector{Float64},
        Prob0::Vector{Float64},
        gamma::Float64
    )
    
    # for e in 1:num_reps
    @inbounds @fastmath for e in 1:num_reps

        #display("e = $e")

        # 1) Find node vjmax with largest D
        vjmax = find_vjmax(Dn)
        if vjmax == 0
            rbp_not_converged = false
            break # i.e., BP has converged
        end        
            
        if Dn[vjmax] < gamma            
            cimax, vjmax = findmaxedge(Residues,alpha,Nc)
        else
            cimax = find_cimax(Residues,Nv,vjmax)
        end
        if cimax == 0.0
            rbp_not_converged = false
            break # i.e., BP has converged
        end

        #2), 3), 4) C2V update
        Nvjmax = Nv[vjmax]
        for ci in Nvjmax
            RBP_update_Lr!(Lr,newLr,Lq,Residues,false,nothing,0.0,Nc,ci,vjmax,msum2)
        end
        Ld = calc_Ld(vjmax,Nvjmax,Lf,Lr)
        bitvector[vjmax] = signbit(Ld)

        Prob0[vjmax] = calc_prob(Ld)
        Dn[vjmax] = 0.0

        maxalp = 0.0
        for vj in Nc[cimax]
            residue = Residues[cimax,vj]
            if residue > maxalp
                maxalp = residue
            end
        end
        alpha[cimax] = maxalp

        for ci in Nvjmax
            # 5) update Nv messages Lq[ci,vnmax]
            li = LinearIndices(Lq)[ci,vjmax]
            alp = Residues[li]
            Lq[li] = tanhLq(Ld,Lr[li],msum_factor)
            # 6) calculate Residues
            Nci = Nc[ci]
            for vj in Nci
                if vj ≠ vjmax
                    li = LinearIndices(Lr)[ci,vj]
                    newlr = calc_Lr(Nci,ci,vj,Lq,msum_factor)
                    newLr[li] = newlr
                    residue = abs(newlr - Lr[li])
                    if residue > alp
                        alp = residue
                    end
                    Residues[li] = residue
                    ln = calc_Ld(vj,Nv[vj],Lf,newLr)
                    new_prob0 = calc_prob(ln)
                    Dn[vj] = abs(Prob0[vj] - new_prob0)                        
                end
            end
            alpha[ci] = alp
        end
    end

    return rbp_not_converged
end

function 
    RBP_update_Lr!(
        Lr::Matrix{Float64},
        newLr::Matrix{Float64},
        Lq::Matrix{Float64},
        Residues::Matrix{Float64},
        RBP_decay::Bool,
        Factors::Union{Matrix{Float64},Nothing},
        decayfactor::Float64,
        Nc::Vector{Vector{Int}},
        ci::Int,
        vjmax::Int,
        msum2::Bool
    )

    # begin
    @fastmath @inbounds begin
        li = LinearIndices(Lr)[ci,vjmax]
        # 2) Decay the RBP factor corresponding to the maximum residue
        if RBP_decay
            Factors[li] *= decayfactor
        end
        # 3) update check to node message Lr[ci,vjmax]
        if msum2
            Lr[li] = calc_Lr_no_opt(Nc[ci],ci,vjmax,Lq)
        else
            Lr[li] = newLr[li]
        end
        # 4) set maximum residue to zero
        Residues[li] = 0.0
    end
end

function 
    calc_Dn!(
        Dn::Vector{Float64},
        Prob0::Vector{Float64},
        newLr::Matrix{Float64},
        Lf::Vector{Float64},
        Nv::Vector{Vector{Int}}
    )

    @inbounds @fastmath for vj in eachindex(Nv)
        ln = calc_Ld(vj,Nv[vj],Lf,newLr)
        new_prob0 = calc_prob(ln)
        Dn[vj] = abs(Prob0[vj] - new_prob0)
    end

end

function 
    init_Prob0!(
       Prob0::Vector{Float64}, 
       Lf::Vector{Float64}
    )

    @inbounds for vj in eachindex(Lf)
        Prob0[vj] = calc_prob(Lf[vj])
    end

end

function calc_prob(ln::Float64)
    @fastmath begin
        exp_ln = exp(ln)
        return exp_ln/(1 + exp_ln)
    end
end

function find_vjmax(Dn::Vector{Float64})

    max_dn = 0.0
    vjmax = 0
    @inbounds @fastmath for vj in eachindex(Dn)
        dn = Dn[vj]
        if dn > max_dn
            max_dn = dn
            vjmax = vj
        end
    end

    return vjmax

end

function 
    find_cimax(
        Residues::Matrix{Float64},
        Nv::Vector{Vector{Int}},
        vjmax::Int
    )

    maxresidue = 0.0
    cimax = 0
    @inbounds @fastmath for ci in Nv[vjmax]
        residue = Residues[ci,vjmax]
        if residue > maxresidue
            maxresidue = residue
            cimax = ci
        end
    end

    return cimax

end


