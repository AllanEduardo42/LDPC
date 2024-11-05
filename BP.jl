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

function 
    BP!(
        supermode::String,
        mode::String,
        stop::Bool,
        test::Bool,
        maxiter::Integer,
        syndrome::Vector{Bool},
        bitvector::Vector{Bool},
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
        Residues::Union{Matrix{<:AbstractFloat},Nothing},
        maxresidue::AbstractFloat,
        maxcoords::Union{Vector{<:Integer},Nothing},
        Factors::Union{Matrix{<:AbstractFloat},Nothing},
        rbpfactor::Union{AbstractFloat,Nothing},
        num_edges::Union{Integer,Nothing},
        Ldn::Union{Vector{<:AbstractFloat},Nothing},
        visited_vns::Union{Vector{Bool},Nothing},
        samples::Union{Vector{<:Integer},Nothing},
        rgn_sample::Union{AbstractRNG,Nothing},
        listsize::Integer,
        listres::Union{Vector{<:AbstractFloat},Nothing},
        listadd::Union{Matrix{<:Integer},Nothing},
        inlist::Union{Matrix{<:Integer},Nothing}
    )
    
    for i in 1:maxiter

        ilbp = (mode == "iLBP" && i > 1)

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
                 ilbp)   
        elseif supermode == "RBP"
            RBP!(Residues,
                 bitvector,
                 Lr,
                 Ms,
                 maxresidue,
                 maxcoords,
                 Lq,
                 Lf,
                 cn2vn,
                 vn2cn,
                 signs,
                 Factors,
                 rbpfactor,
                 num_edges,
                 Ldn,
                 samples,
                 rgn_sample,
                 listsize,
                 listres,
                 listadd,
                 inlist)
            # reset factors
            resetfactors!(Factors,vn2cn)
        end

        ilbp ? nothing : calc_syndrome!(syndrome,bitvector,cn2vn)

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
                elseif ilbp
                    syndrome .= true
                end
                if stop
                    if i < maxiter
                        @inbounds decoded[i+1:end] .= true
                    end
                    break
                end
            end
            biterror .= (bitvector .≠ codeword)
            @inbounds ber[i] = sum(biterror)
        end
    end

end