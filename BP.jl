################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Main loop of the BP routine

include("update_Lr.jl")
include("update_Lq.jl")
include("minsum.jl")
include("flooding.jl")
include("calc_syndrome.jl")
include("LBP.jl")
# include("LBP_NO_OPT.jl")
include("RBP.jl")
include("local_RBP.jl")

function 
    BP!(
        address::Union{Matrix{<:Integer},Nothing},
        addressinv::Union{Matrix{<:Integer},Nothing},
        supermode::String,
        stop::Bool,
        test::Bool,
        maxiter::Integer,
        syndrome::Vector{Bool},
        syndrome2::Vector{Bool},
        bitvector::Vector{Bool},
        bitvector2::Vector{Bool},
        codeword::Vector{Bool},
        biterror::Vector{Bool},
        ber::Vector{<:Integer},
        decoded::Vector{Bool},
        Lf::Array{<:AbstractFloat},
        Lq::Array{<:AbstractFloat},
        Lr::Array{<:AbstractFloat},
        Ms::Union{Matrix{<:AbstractFloat},Nothing},      
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        Lrn::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},        
        phi::Union{Vector{<:AbstractFloat},Nothing},
        printtest::Union{Bool,Nothing},        
        residues::Union{Vector{<:AbstractFloat},Nothing},
        Factors::Union{Matrix{<:AbstractFloat},Nothing},
        decay::Union{AbstractFloat,Nothing},
        num_edges::Union{Integer,Nothing},
        Ldn::Union{Vector{<:AbstractFloat},Nothing},
        visited_vns::Union{Vector{Bool},Nothing},
        rgn_sample::Union{AbstractRNG,Nothing},
        listsize::Integer,
        listsize2::Integer,
        listres1::Union{Vector{<:AbstractFloat},Nothing},
        listm1::Union{Vector{<:Integer},Nothing},
        listn1::Union{Vector{<:Integer},Nothing},
        listres2::Union{Vector{<:AbstractFloat},Nothing},
        listm2::Union{Vector{<:Integer},Nothing},
        listn2::Union{Vector{<:Integer},Nothing},
        inlist::Union{Matrix{<:Integer},Nothing},
        maxresidues::Vector{<:AbstractFloat},
        maxcoords::Vector{<:Integer},
        maxcoords_alt::Vector{<:Integer}
    )
    
    for i in 1:maxiter

        if test && printtest  
            println("### Iteration #$i ###")
        end

        if supermode == "Flooding"
            flooding!(bitvector,
                      Lq,
                      Lr,
                      Lf,
                      cn2vn,
                      vn2cn,
                      Lrn,
                      signs,
                      phi)  
        elseif supermode == "LBP"
            LBP!(bitvector,
                 Lq,
                 Lr,
                 Lf,
                 cn2vn,
                 vn2cn,
                 Lrn,
                 syndrome,
                 Ldn,
                 visited_vns,
                 i)   
        elseif supermode == "RBP"
            RBP!(address,
                 addressinv,
                 residues,
                 bitvector,
                 bitvector2,
                 Lr,
                 Ms,
                 Lq,
                 Lf,
                 cn2vn,
                 vn2cn,
                 Lrn,
                 signs,
                 phi,
                 Factors,
                 decay,
                 num_edges,
                 Ldn,
                 rgn_sample,
                 listsize,
                 listsize2,
                 listres1,
                 listm1,
                 listn1,
                 listres2,
                 listm2,
                 listn2,
                 inlist,
                 syndrome,
                 syndrome2)
            # reset factors
            resetfactors!(Factors,vn2cn)
        elseif supermode == "Local-RBP"
            local_RBP!(
                maxresidues,
                maxcoords,
                bitvector,
                Lr,
                Ms,
                Lq,
                Lf,
                cn2vn,
                vn2cn,
                Lrn,
                signs,
                phi,
                Factors,
                decay,
                num_edges,
                Ldn,
                rgn_sample,
                maxcoords_alt
            )
            # reset factors
            resetfactors!(Factors,vn2cn)
        end

        calc_syndrome!(syndrome,bitvector,cn2vn)

        if test && printtest    
                println("Max LLR estimate errors: ")
                for j in eachindex(bitvector)
                    print(Int(bitvector[j] != codeword[j]))
                    if j%80 == 0
                        println()
                    end
                end     
                println() 
                println("Syndrome: ")
                for j in eachindex(syndrome)
                    print(Int(syndrome[j]))
                    if j%80 == 0
                        println()
                    end
                end     
                println()
                if iszero(syndrome) && stop
                    break
                end
                println()     
        else
            if iszero(syndrome)
                if bitvector == codeword
                    @inbounds decoded[i] = true
                end
                if stop
                    if i < maxiter
                        @inbounds decoded[i+1:end] .= decoded[i]
                    end
                    break
                end
            end
            biterror .= (bitvector .≠ codeword)
            @inbounds ber[i] = sum(biterror)
        end
    end

end