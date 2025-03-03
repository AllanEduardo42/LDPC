################################################################################
# Allan Eduardo Feitosa
# 3 Mar 2025
# RBP Sum-Product Algorithm using the only local strategy with residual decaying
# factor

include("./RBP functions/RBP_update_Lr.jl")
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
        Ms::Matrix{<:AbstractFloat},
        Factors::Matrix{<:AbstractFloat},        
        all_max_res_alt::Union{Vector{<:AbstractFloat},Nothing},
        test::Bool,
        rng_rbp::AbstractRNG,
        max_residue::Vector{<:AbstractFloat},
        max_indices::Vector{<:Integer}
    )

    @fastmath @inbounds for e in 1:num_edges

        # display("e = $e")
        # display(max_residue)

        # if in test mode, store the values of the maximum residues
        if test
            all_max_res_alt[e] = max_residue[1] 
        end 

        # 1) get the largest residues coordenates
        if max_residue[1] == 0.0
            cnmax = rand(rng_rbp,1:length(cn2vn))
            vnmax = rand(rng_rbp,cn2vn[cnmax])
            lmax = LinearIndices(Factors)[cnmax,vnmax]
        else
            lmax = max_indices[1]
            ci = CartesianIndices(Factors)[lmax]
            cnmax, vnmax = ci[1], ci[2]
        end        

        # 2) Decay the RBP factor corresponding to the maximum residue
        Factors[lmax] *= decayfactor

        # 3) update check to node message Lr[cnmax,vnmax]
        RBP_update_Lr!(lmax,Lr,Ms,cnmax,vnmax,cn2vn,Lq,Lrn,signs,phi)
        
        # 4) clear the max residue
        max_residue[1] = 0.0
        # max_residue[2] = 0.0

        # 5) update vn2cn messages Lq[vnmax,m] and bitvector[vnmax]
        bitvector[vnmax] = update_Lq!(Lq,Lr,Lf[vnmax],vnmax,vn2cn,Lrn)

        # 6) calculate residues
        for m in vn2cn[vnmax]
            if m ≠ cnmax    
                # calculate the new check to node messages
                update_Lr!(Ms,Lq,m,cn2vn,Lrn,signs,phi)
                # calculate the residues
                for n in cn2vn[m]
                    if n ≠ vnmax
                        residue, li = calc_residue(Ms,Lr,Factors,Lrn,Lq,m,n)
                        if residue > max_residue[1]
                            max_residue[2] = max_residue[1]
                            max_residue[1] = residue
                            max_indices[2] = max_indices[1]
                            max_indices[1] = li
                        elseif residue > max_residue[2]
                            max_residue[2] = residue
                            max_indices[2] = li
                        end   
                    end
                end
            end
        end
    
        # update list
        if max_residue[1] < max_residue[3]
            max_residue[1], max_residue[3] = max_residue[3], max_residue[1] 
            max_indices[1], max_indices[3] = max_indices[3], max_indices[1]    
        else
            max_residue[3] = max_residue[2]
            max_indices[3] = max_indices[2]
        end
    end
end

