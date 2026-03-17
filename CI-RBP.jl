################################################################################
# Allan Eduardo Feitosa
# 16 Jan 2026
# CI-RBP Algorithm

include("./RBP functions/findmaxedge.jl")

function
    CI_RBP!(
        bitvector::Vector{Bool},
        V2C::Matrix{Float64},
        C2V::Matrix{Float64},
        prior_LLRs::Vector{Float64},
        Nc::Vector{Vector{Int}},
        Nv::Vector{Vector{Int}},
        phi::Union{Vector{Float64},Nothing},
        msum_factor::Union{Float64,Nothing},
        msum2::Bool,
        num_reps::Int,
        newC2V::Matrix{Float64},
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

        limax = LinearIndices(C2V)[cimax,vjmax]
        # 2) update check to node message C2V[ci,vjmax]
        if msum2
            C2V[limax] = calc_C2V_no_opt(Nc[cimax],cimax,vjmax,V2C)
        else
            C2V[limax] = newC2V[limax]
        end
        # 3) set maximum residue to zero
        Residues[limax] = 0.0

        # 4) update LLR[vjmax]
        Nvjmax = Nv[vjmax]
        post_LLR = calc_post_LLR(vjmax,Nvjmax,prior_LLRs,C2V)
        bitvector[vjmax] = signbit(post_LLR)

        Prob0[vjmax] = calc_prob(post_LLR)
        Dn[vjmax] = 0.0

        # update alpha[cimax]
        maxalp = 0.0
        for vj in Nc[cimax]
            residue = Residues[cimax,vj]
            if residue > maxalp
                maxalp = residue
            end
        end
        alpha[cimax] = maxalp

        for ci in Nvjmax
            if ci ≠ cimax
                # 5) update Nv messages V2C[ci,vnmax]
                li = LinearIndices(V2C)[ci,vjmax]
                alp = Residues[li]
                V2C[li] = tanh_V2C(post_LLR,C2V[li],msum_factor)
                # 6) calculate Residues
                Nci = Nc[ci]
                for vj in Nci
                    if vj ≠ vjmax
                        li = LinearIndices(C2V)[ci,vj]
                        newc2v = calc_C2V(Nci,ci,vj,V2C,msum_factor)
                        newC2V[li] = newc2v
                        residue = abs(newc2v - C2V[li])
                        if residue > alp
                            alp = residue
                        end
                        Residues[li] = residue
                        calc_Dn!(Dn,Prob0,newC2V,prior_LLRs,vj,Nv)                    
                    end
                end
                alpha[ci] = alp
            end
        end
    end

    return rbp_not_converged
end

function 
    calc_Dn!(
        Dn::Vector{Float64},
        Prob0::Vector{Float64},
        newC2V::Matrix{Float64},
        prior_LLRs::Vector{Float64},
        vj::Int,
        Nv::Vector{Vector{Int}}
    )

    @inbounds @fastmath begin
        llr = calc_post_LLR(vj,Nv[vj],prior_LLRs,newC2V)
        new_prob0 = calc_prob(llr)
        Dn[vj] = abs(Prob0[vj] - new_prob0)
    end


end

function calc_prob(llr::Float64)
    @fastmath begin
        exp_llr = exp(llr)
        return exp_llr/(1 + exp_llr)
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


