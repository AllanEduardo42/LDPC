################################################################################
# Allan Eduardo Feitosa
# 11 out 2024
# RBP Sum-Product Algorithm

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/RBP_set_zero_or_remove.jl")
include("./RBP functions/calc_residues.jl")
include("./RBP functions/findmaxcoords.jl")
include("./RBP functions/update_list.jl")
include("calc_syndrome.jl")

#RBP
function
    RBP!(
        H::BitMatrix,
        address::Union{Matrix{<:Integer},Nothing},
        addressinv::Union{Matrix{<:Integer},Nothing},
        residues::Union{Vector{<:AbstractFloat},Nothing},
        bitvector::Vector{Bool},
        bitvector2::Vector{Bool},
        Lr::Matrix{<:AbstractFloat},
        Ms::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        Lrn::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},        
        phi::Union{Vector{<:AbstractFloat},Nothing},
        Factors::Matrix{<:AbstractFloat},
        decayfactor::AbstractFloat,
        num_edges::Integer,
        Ldn::Vector{<:AbstractFloat},
        rng_sample::AbstractRNG,
        listsize::Integer,
        listsize2::Integer,
        listres::Vector{<:AbstractFloat},
        listm::Vector{<:Integer},
        listn::Vector{<:Integer},
        listres2::Union{Vector{<:AbstractFloat},Nothing},
        listm2::Union{Vector{<:Integer},Nothing},
        listn2::Union{Vector{<:Integer},Nothing},
        inlist::Union{Matrix{<:Integer},Nothing},
        syndrome::Vector{Bool},
        syndrome2::Vector{Bool}
    )

    # sum_syndrome = zeros(Int,4)
    # sum_syndrome2 = zeros(Int,listsize)

    @fastmath for e in 1:num_edges

        # println("e = $e")
        # x = e +(iter-1)*num_edges
        # global Residues[x] = mean(residues)

        # display([listm listn listres])

        # sum_syndrome .= M

        # calc_syndrome!(syndrome,bitvector,cn2vn)
        # # println("sum syndrome = $(sum(syndrome))")
        # sum_syndrome .= 0

        # for i=1:4
        #     if listres[i] != 0
        #         cnmax = listm[i]
        #         vnmax = listn[i]
        #     else
        #         break
        #     end
        #     bitvector2 = copy(bitvector)
        #     # syndrome2 = copy(syndrome)
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
        #     # bitvector2[vnmax] = new_bit           
        #     if bitvector[vnmax] != new_bit
        #         for m in eachindex(cn2vn)
        #             if H[nl+m]
        #                 if syndrome[m] == true
        #                     sum_syndrome[i] -= 1
        #                 else
        #                     sum_syndrome[i] += 1
        #                 end
        #             end
        #         end
        #     end 
        #     # calc_syndrome!(syndrome2,bitvector2,cn2vn)     
        #     # sum_syndrome2[i] = sum(syndrome2)     
        # end

        # # println("[index sum_syndrome listres] =")
        # # display([collect(1:4) sum_syndrome listres[1:4]])

        # _,index = findmin(sum_syndrome)

        # display(index)

        #######

        # cum = zeros(listsize+1)

        # for i=1:listsize
        #     if listres[i] != 0
        #         cnmax = listm[i]
        #         vnmax = listn[i]
        #     else
        #         break
        #     end
        #     cum[i] = Lf[vnmax]
        #     nl = LinearIndices(Lr)[1,vnmax]-1
        #     for m in vn2cn[vnmax]
        #         if m == cnmax
        #             cum[i] += Ms[nl+m]
        #         else
        #             cum[i] += Lr[nl+m]
        #         end
        #     end
        # end

        # # display(cum)

        # _,index = findmax(abs.(cum))

        # display(index)

        index = 1
        if listres[index] != 0
            cnmax = listm[index]
            vnmax = listn[index]
        else
            cnmax = rand(rng_sample,1:length(cn2vn))
            vnmax = rand(rng_sample,cn2vn[cnmax])
        end

        # display([cnmax vnmax])

        # display([listm listn listres cum])

        # 2) for optimization, get the linear indice of the maximum residue
        lmax = LinearIndices(Factors)[cnmax,vnmax]

        # 3) Decay the RBP factor corresponding to the maximum residue
        Factors[lmax] *= decayfactor

        # 4) update check to node message Lr[cnmax,vnmax]
        RBP_update_Lr!(lmax,cnmax,vnmax,cn2vn,Lq,Lr,Ms,Lrn,signs,phi)

        # 5) set maximum residue to zero or remove it from the list
        set_zero_or_remove!(addressinv,residues,lmax,listsize,listres,listm,
                            listn,inlist,index)

        # 6) update Ldn[vmax] and bitvector[vnmax]
        # Ldn[vnmax] = Lf[vnmax]
        # nl = LinearIndices(Lr)[1,vnmax]-1
        # for m in vn2cn[vnmax]
        #     Ldn[vnmax] += Lr[nl+m]
        # end
        # bitvector[vnmax] = signbit(Ldn[vnmax])
        _, bitvector[vnmax] = update_Lq!(Lq,Lr,Lf[vnmax],vnmax,vn2cn,Lrn)

        # 7) update vn2cn messages Lq[vnmax,m], ∀m ≠ cnmax, and calculate residues
        leaf = true # suppose node vnmax is a leaf in the graph
        count_size = listsize2
        # display(vn2cn[vnmax])
        for m in vn2cn[vnmax]
            if m ≠ cnmax
                leaf = false # vnmax is not a leaf
                # Lq[vnmax,m] = Ldn[vnmax] - Lr[nl+m]
                count_size = calc_residues!(addressinv,residues,Factors,Ms,Lr,Lq,Lrn,signs,
                               phi,vnmax,m,cn2vn,listres,listm,listn,listres2,
                               listm2,listn2,listsize,listsize2,count_size,inlist)
            end
        end

        # 8) if vnmax is a leaf in the graph
        if leaf
            Lq[vnmax,cnmax] = Ldn[vnmax] - Lr[nl+cnmax]
            find_local_maxresidue!(maxresidues,Factors,Ms,Lr,Lq,
                Lrn,signs,phi,vnmax,cnmax,cn2vn,maxcoords)
        end

        # display([listm2 listn2 listres2])

        # 9) Update list 1
        if listsize2 ≠ 0
            count = 0
            k = listsize
            while count < listsize2-1
                count += 1
                k -= 1
                m = listm[k]
                if m != 0
                    listres[k] = 0.0
                    n = listn[k]
                    inlist[m,n] = false
                end
            end
            # display([listm listn listres])
            for k in count_size:-1:1
                update_list!(inlist,listres,listm,listn,listres2[k],listm2[k],
                             listn2[k],listsize)
            end
            listres2 .*= 0.0
            listm2 .*= 0
            listn2 .*= 0
        end

        # display([listm listn listres])

        # 10) find maximum residue (only for the original RBP)

        findmaxcoords!(address,residues,listres,listm,listn)

    end
end

