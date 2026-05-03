################################################################################
# Allan Eduardo Feitosa
# 07 Apr 2025
# RBP and List-RBP Algorithm with residual decaying factor

include("./List functions/update_main_list.jl")
include("./List functions/add_to_list.jl")
include("./List functions/remove_from_list.jl")
include("./List functions/find_list_pos.jl")

function
    List_RBP!(
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
        Residuals::Matrix{Float64},
        decayfactor::Float64,
        Factors::Matrix{Float64},
        coords::Matrix{Int},
        inlist::Matrix{Bool},        
        list::Vector{Float64},
        local_list::Union{Vector{Float64},Nothing},
        local_coords::Union{Matrix{Int},Nothing},
        listsize::Int,
        listsize2::Int        
    )

    rbp_not_converged = true

    @fastmath @inbounds for e in 1:num_reps

    # display("e = $e")

       # 1) Find largest residual and coordenates

       # if max residual is equal to zero, refill the list
        if list[1] == 0.0
            for ci in eachindex(Nc)
                for vj in Nc[ci]
                    li = LinearIndices(Residuals)[ci,vj]
                    residual = Residuals[li]
                    # add residual to main list
                    add_to_list!(inlist,list,coords,residual,li,ci,vj,listsize)
                end
            end
            # if max residual is still zero, List-RBP has converged
            if list[1] == 0.0
                rbp_not_converged = false
                break
            end
        end
        cimax = coords[1,1]
        vjmax = coords[2,1]
        limax = coords[3,1]

        # 2) update C2V message C2V[cimax,vjmax]
        C2V[limax] = newC2V[limax]

        # 3) set maximum residual to zero
        Residuals[limax] = 0.0
        remove_from_list!(limax,listsize,list,coords,inlist,1)

        # 4) Decay the RBP factor corresponding to the maximum residual
        Factors[limax] *= decayfactor

        Nvjmax = Nv[vjmax]
        post_LLR = calc_post_LLR(vjmax,Nvjmax,prior_LLRs,C2V)
        bitvector[vjmax] = signbit(post_LLR)

        for ci in Nvjmax
            if ci ≠ cimax
                # 5) update Nv messages V2C[ci,vnmax]
                li = LinearIndices(V2C)[ci,vjmax]
                V2C[li] = tanh_V2C(post_LLR,C2V[li],msum_factor)
                # 6) calculate residues
                Nci = Nc[ci]
                for vj in Nci
                    if vj ≠ vjmax
                        li = LinearIndices(C2V)[ci,vj]
                        newc2v = calc_C2V(Nci,ci,vj,V2C,msum_factor)
                        newC2V[li] = newc2v
                        residual = abs(newc2v - C2V[li])
                        residual *= Factors[li]
                        Residuals[li] = residual
                        # if egde "li" is on the main list, find position and remove
                        if inlist[li]
                            pos = find_list_pos(li,listsize,coords,ci,vj)
                            remove_from_list!(li,listsize,list,coords,inlist,pos)
                        end
                        # add residual to local list
                        add_to_list!(nothing,local_list,local_coords,residual,li,
                                                                ci,vj,listsize2)
                    end
                end
            end
        end
        update_main_list!(list,coords,local_list,local_coords,listsize,listsize2,
                                                                         inlist)
    end

    return rbp_not_converged
end

