################################################################################
# Allan Eduardo Feitosa
# 31 Mar 2025
# RBP Sum-Product Algorithm with residual decaying factor

include("./RBP functions/findmaxnode.jl")
include("./RBP functions/calc_residue.jl")

function
    LD_RBP!(
        bitvector::Vector{Bool},
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},
        Lf::Vector{Float64},
        Nc::Vector{Vector{Int}},
        Nv::Vector{Vector{Int}},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{Float64},Nothing},
        decayfactor::Float64,
        num_reps::Int,
        newLr::Matrix{Float64},
        Factors::Vector{Float64},
        alpha::Vector{Float64},
        rbp_not_converged::Bool
    )

    @fastmath @inbounds for e in 1:num_reps

        # display("e = $e")

        vjmax = findmaxnode(alpha)
        if vjmax == 0
            rbp_not_converged = false
            break # i.e., BP has converged
        end
        Factors[vjmax] *= decayfactor
        alpha[vjmax] = 0.0

        Nvjmax = Nv[vjmax]
        for ci in Nvjmax
            li = LinearIndices(Lr)[ci,vjmax]
            Lr[li] = newLr[li]
        end 

        Ld = calc_Ld(vjmax,Nvjmax,Lf,Lr)
        bitvector[vjmax] = signbit(Ld)

        for ci in Nvjmax
            # 5) update Nv messages Lq[ci,vnmax]
            li = LinearIndices(Lq)[ci,vjmax]
            Lq[li] = tanh(0.5*(Ld - Lr[li]))
            # 6) calculate alpha
            Nci = Nc[ci]
            for vj in Nci
                if vj ≠ vjmax
                    li = LinearIndices(Lr)[ci,vj]
                    newLr[li] = calc_Lr(Nci,ci,vj,Lq)
                    alp = 0.0
                    for ca in Nv[vj]
                        li = LinearIndices(Lr)[ca,vj]
                        alp += newLr[li] - Lr[li]
                    end
                    alpha[vj] = abs(alp)
                end
            end
        end
    end  
    
    return rbp_not_converged
    
end