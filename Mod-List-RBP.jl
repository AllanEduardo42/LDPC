################################################################################
# Allan Eduardo Feitosa
# 24 Feb 2024
# Modified List-RBP Sum-Product Algorithm with residual decaying factor

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/calc_residues.jl")
include("./RBP functions/update_list.jl")
include("./RBP functions/decay.jl")
include("./RBP functions/remove_from_list.jl")

function
    mod_list_RBP!(
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
        Factors::Matrix{<:AbstractFloat},        
        all_max_res_alt::Union{Vector{<:AbstractFloat},Nothing},
        test::Bool,
        rng_rbp::AbstractRNG,
        listsizes::Vector{<:Integer},
        listres::Vector{<:AbstractFloat},
        listm::Vector{<:Integer},
        listn::Vector{<:Integer},
        listres2::Union{Vector{<:AbstractFloat},Nothing},
        listm2::Union{Vector{<:Integer},Nothing},
        listn2::Union{Vector{<:Integer},Nothing},
        inlist::Union{Matrix{<:Integer},Nothing},
        syndrome::Vector{Bool}        
    )

    # @inbounds sum_syndrome = zeros(Int,listsizes[3])

    @inbounds @fastmath for e in 1:2

        display("e = $e")
        # calc_syndrome!(syndrome,bitvector,cn2vn)
        # Syndrome[e] = sum(syndrome)
        # display("sum syndrome = $(sum(syndrome))")

        if test
            all_max_res_alt[e] = listres[1]     
        end

        # sum_syndrome .= 0
        # for i=1:listsizes[3]
        #     if listres[i] != 0
        #         cnmax = listm[i]
        #         vnmax = listn[i]
        #     else
        #         break
        #     end
        #     nl = LinearIndices(Lr)[1,vnmax]-1
        #     cum = Lf[vnmax]
        #     for m in vn2cn[vnmax]
        #         if m == cnmax
        #             cum += Ms[nl+m]
        #         else
        #             cum += Lr[nl+m]
        #         end
        #     end
        #     new_bit = signbit(cum)         
        #     if bitvector[vnmax] != new_bit
        #         for m in vn2cn[vnmax]
        #             if syndrome[m] == true
        #                 sum_syndrome[i] -= 1
        #             else
        #                 sum_syndrome[i] += 1
        #             end
        #         end
        #     end
        #     # bitvector2 = copy(bitvector)
        #     # bitvector2[vnmax] = new_bit
        #     # calc_syndrome!(syndrome,bitvector2,cn2vn)
        #     # sum_syndrome[i] = sum(syndrome)
        # end

        # display(sum_syndrome)

        # minsyn, index = findmin(sum_syndrome)

        # display("index = $index")

        cum = zeros(listsizes[3])
        diff_cum = zeros(listsizes[3])
        for i=1:listsizes[3]
            if listres[i] != 0
                cnmax = listm[i]
                vnmax = listn[i]
            else
                break
            end
            cum[i] = Lf[vnmax]
            nl = LinearIndices(Lr)[1,vnmax]-1
            for m in vn2cn[vnmax]
                if m == cnmax
                    cum[i] += Ms[nl+m]
                else
                    cum[i] += Lr[nl+m]
                end
            end
            diff_cum[i] = abs(cum[i] - (Lq[vnmax,cnmax]+Lr[cnmax,vnmax]))
        end

        display(diff_cum)
        display(listres)

        _,index = findmax(abs.(diff_cum))

        display(index)
        
        # 1) get the largest residues coordenates if not clipped
        new_listsize2 = listsizes[2]
        # index = 1
        # if listsizes[1] > 1
        #     search_index = true
        #     i = 0
        #     while search_index
        #         i += 1
        #         if abs(listres[index]) < CLIP
        #             search_index = false
        #             index = i
        #         else
        #             if listsizes[2] == 1
        #                 new_listsize2 += 1
        #             end
        #             m = listm[index]
        #             n = listn[index]
        #             lmax = LinearIndices(Factors)[m,n]
        #             remove_from_list!(lmax,listsizes[1],listres,listm,listn,inlist,index)
        #         end
        #     end
        # end
        if listres[index] != 0
            cnmax = listm[index]
            vnmax = listn[index]
        else
            cnmax = rand(rng_rbp,1:length(cn2vn))
            vnmax = rand(rng_rbp,cn2vn[cnmax])
        end        

        # 2) Decay the RBP factor corresponding to the maximum residue
        lmax = decay!(cnmax,vnmax,Factors,decayfactor)

        # 3) update check to node message Lr[cnmax,vnmax]
        RBP_update_Lr!(lmax,Lr,Ms,cnmax,vnmax,cn2vn,Lq,Lrn,signs,phi)

        # 4) Remove max residue from the list and update the list
        remove_from_list!(lmax,listsizes[1],listres,listm,listn,inlist,index) 

        # 5) update Ldn[vmax] and bitvector[vnmax]
        bitvector[vnmax] = update_Lq!(Lq,Lr,Lf[vnmax],vnmax,vn2cn,Lrn)
        # Obs: this code update Lq[vnmax,cnmax]: no difference in performance 
        #      was detected. The original implementation, which doesn't update
        #      Lq[vnmax,cnmax], is as follows (steps 1, 2 and 3):
        # step 1) Ldn[vnmax], nl = calc_Ld(vnmax,vn2cn,Lf[vnmax],Lr)
        # step 2) bitvector[vnmax] = signbit(Ldn[vnmax]
        # step 3) see below.
        # calc_syndrome!(syndrome, bitvector,cn2vn)
        # display("new sum syndrome = $(sum(syndrome))")

        # 6) update vn2cn messages Lq[vnmax,m], ∀m ≠ cnmax, and calculate residues
        checks = vn2cn[vnmax]
        if length(checks) > 1
            for m in checks
                if m ≠ cnmax
                    # step 3) Lq[vnmax,m] = Ldn[vnmax] - Lr[nl+m]        
                    new_listsize2 = calc_residues!(Factors,Ms,Lr,Lq,Lrn,signs,
                            phi,vnmax,m,cn2vn,listres,listm,listn,listres2,
                            listm2,listn2,listsizes,new_listsize2,inlist)
                end
            end
        else # if vnmax is a leaf in the graph
            new_listsize2 = calc_residues!(Factors,Ms,Lr,Lq,Lrn,signs,phi,vnmax,
                            m,cn2vn,listres,listm,listn,listres2,listm2,listn2,
                            listsizes,new_listsize2,inlist)
        end

        # 7) update list
        if listsizes[2] ≠ 0
            count = 0
            k = listsizes[1]
            while count < listsizes[2]-1
                count += 1
                k -= 1
                m = listm[k]
                if m != 0
                    listres[k] = 0.0
                    n = listn[k]
                    inlist[m,n] = false
                end
            end
            for k in new_listsize2:-1:1
                update_list!(inlist,listres,listm,listn,listres2[k],listm2[k],
                             listn2[k],listsizes[1])
            end
            listres2 .*= 0.0
            listm2 .*= 0
            listn2 .*= 0
        end     

    end
end

