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

    sum_syndrome = zeros(Int,listsize)

    @fastmath @inbounds for e in 1:num_edges

        # println("e = $e")

        # sum_syndrome .= M

        # for i=1:listsize
        #     if listres[i] != 0
        #         vnmax = listn[i]
        #     else
        #         break
        #     end
        #     bitvector2 = copy(bitvector)
        #     syndrome2 = copy(syndrome)
        #     nl = LinearIndices(Lr)[1,vnmax]-1
        #     cum = Lf[vnmax]
        #     for m in vn2cn[vnmax]
        #         cum += Ms[nl+m]
        #     end
        #     bitvector2[vnmax] = signbit(cum)
        #     for m in eachindex(cn2vn)
        #         for n in cn2vn[m]
        #             if n == vnmax
        #                 syndrome2[m] = syndrome[m] ⊻ bitvector[vnmax] ⊻ bitvector2[vnmax]
        #             else
        #                 syndrome2[m] = syndrome[m]
        #             end
        #         end
        #     end
        #     # calc_syndrome!(syndrome2,bitvector2,cn2vn)
        #     sum_syndrome[i] = sum(syndrome2)
        # end

        # display([sum_syndrome listres[1:end-1]])

        # _,index = findmin(sum_syndrome)

        # display(index)

        # 1) verify if the list was not updated
        index = 1
        if listres[index] != 0
            cnmax = listm[index]
            vnmax = listn[index]
        else
            cnmax = rand(rng_sample,1:length(cn2vn))
            vnmax = rand(rng_sample,cn2vn[cnmax])
        end

        # 2) for optimization, get the linear indice of the maximum residue
        lmax = LinearIndices(Factors)[cnmax,vnmax]

        # 3) Decay the RBP factor corresponding to the maximum residue
        Factors[lmax] *= decayfactor

        # 4) update check to node message Lr[cnmax,vnmax]
        RBP_update_Lr!(lmax,cnmax,vnmax,cn2vn,Lq,Lr,Ms,Lrn,signs,phi)

        # 5) set maximum residue to zero or remove it from the list
        set_zero_or_remove!(addressinv,residues,lmax,listsize,listres,listm,
                            listn,inlist,1)

        # 6) update Ldn[vmax] and bitvector[vnmax]
        Ldn[vnmax] = Lf[vnmax]
        nl = LinearIndices(Lr)[1,vnmax]-1
        for m in vn2cn[vnmax]
            Ldn[vnmax] += Lr[nl+m]
        end
        bitvector[vnmax] = signbit(Ldn[vnmax])

        # 7) update vn2cn messages Lq[vnmax,m], ∀m ≠ cnmax, and calculate residues
        leaf = true # suppose node vnmax is a leaf in the graph
        for m in vn2cn[vnmax]
            if m ≠ cnmax
                leaf = false # vnmax is not a leaf
                Lq[vnmax,m] = Ldn[vnmax] - Lr[nl+m]
                calc_residues!(addressinv,residues,Factors,Ms,Lr,Lq,Lrn,signs,
                               phi,vnmax,m,cn2vn,listres,listm,listn,listres2,
                               listm2,listn2,listsize,listsize2,inlist)
            end
        end

        # 8) if vnmax is a leaf in the graph
        if leaf
            Lq[vnmax,cnmax] = Ldn[vnmax] - Lr[nl+cnmax]
            find_local_maxresidue!(maxresidues,Factors,Ms,Lr,Lq,
                Lrn,signs,phi,vnmax,cnmax,cn2vn,maxcoords)
        end

        # 9) Update list 1
        if listsize2 ≠ 0
            for k in listsize2:-1:1
                update_list!(inlist,listres,listm,listn,listres2[k],listm2[k],
                             listn2[k],listsize)
            end
            listres2 .*= 0.0
            listm2 .*= 0
            listn2 .*= 0
        end

        # 10) find maximum residue (only for the original RBP)

        findmaxcoords!(address,residues,listres,listm,listn)

    end
end

