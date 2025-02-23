################################################################################
# Allan Eduardo Feitosa
# 13 jan 2025
# RBP Sum-Product Algorithm using the only local strategy with residual decaying
# factor

include("./RBP functions/RBP_update_Lr.jl")
include("./RBP functions/find_local_maxresidue.jl")
include("./RBP functions/decay.jl")


function
    local_RBP!(
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
        Ldn::Vector{<:AbstractFloat},
        Ms::Matrix{<:AbstractFloat},
        Factors::Matrix{<:AbstractFloat},        
        all_max_res_alt::Union{Vector{<:AbstractFloat},Nothing},
        test::Bool,
        rng_rbp::AbstractRNG,
        largest_res::Vector{<:AbstractFloat},
        largestcoords::Vector{<:Integer},
        largestcoords_alt::Vector{<:Integer}
    )

    largest_res_alt = largest_res[2]

    @fastmath @inbounds for e in 1:num_edges

        # display("e = $e")

        # if in test mode, store the values of the maximum residues
        if test
            all_max_res_alt[e] = largest_res[1] 
        end 

        # 1) get the largest residues coordenates
        if largest_res[1] == 0
            cnmax = rand(rng_rbp,1:length(cn2vn))
            vnmax = rand(rng_rbp,cn2vn[cnmax])
        else
            largest_res[1] = 0
            cnmax = largestcoords[1]
            vnmax = largestcoords[2]
        end        

        # 2) Decay the RBP factor corresponding to the maximum residue
        lmax = decay!(cnmax,vnmax,Factors,decayfactor)

        # 3) update check to node message Lr[cnmax,vnmax]
        RBP_update_Lr!(lmax,Lr,Ms,cnmax,vnmax,cn2vn,Lq,Lrn,signs,phi)    

        # 4) update Ldn[vmax] and bitvector[vnmax]
        Ldn[vnmax], nl = calc_Ld(vnmax,vn2cn,Lf[vnmax],Lr)
        bitvector[vnmax] = signbit(Ldn[vnmax])

        # 5) update vn2cn messages Lq[vnmax,m], ∀m ≠ cnmax, and calculate residues
        leaf = true # suppose node vnmax is a leaf in the graph
        for m in vn2cn[vnmax]
            if m ≠ cnmax
                leaf = false # vnmax is not a leaf
                Lq[vnmax,m] = Ldn[vnmax] - Lr[nl+m]
                find_local_maxresidue!(largest_res,Factors,Ms,Lr,Lq,
                    Lrn,signs,phi,vnmax,m,cn2vn,largestcoords)
            end
        end

        # 6) if vnmax is a leaf in the graph
        if leaf
            Lq[vnmax,cnmax] = Ldn[vnmax] - Lr[nl+cnmax]
            find_local_maxresidue!(largest_res,Factors,Ms,Lr,Lq,
                Lrn,signs,phi,vnmax,cnmax,cn2vn,largestcoords)
        end

        # 7) update list
        if largest_res[1] < largest_res_alt
            largestcoords[1], largestcoords_alt[1] = largestcoords_alt[1], largestcoords[1]
            largestcoords[2], largestcoords_alt[2] = largestcoords_alt[2], largestcoords[2]
            largest_res[1], largest_res_alt = largest_res_alt, largest_res[1]            
        else
            largest_res_alt = largest_res[2]
            largestcoords_alt[1] = largestcoords[3]
            largestcoords_alt[2] = largestcoords[4]
        end
    end
end

