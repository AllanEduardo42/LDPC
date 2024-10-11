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
        max::Integer,
        syndrome::Vector{Bool},
        d::Vector{Bool},
        c::Vector{Bool},
        bit_error::Vector{Bool},
        ber::Vector{<:AbstractFloat},
        Lf::Array{<:AbstractFloat},
        Lq::Array{<:AbstractFloat},
        Lr::Array{<:AbstractFloat},
        Ms::Union{Matrix{<:AbstractFloat},Nothing},      
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        Lrn::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},        
        phi::Union{Vector{<:AbstractFloat},Nothing},
        printing::Union{Bool,Nothing},        
        Residues::Union{Matrix{<:AbstractFloat},Nothing},
        maxcoords::Union{Vector{<:Integer},Nothing},
        Factors::Union{Matrix{<:AbstractFloat},Nothing},
        H::BitMatrix,
        rbpfactor::Union{AbstractFloat,Nothing},
        num_edges::Union{Integer,Nothing},
        Ldn::Union{Vector{<:AbstractFloat},Nothing},
        visited_vns::Union{Vector{Bool},Nothing},
        samples::Union{Vector{<:Integer},Nothing},
        rgn_sample::Union{AbstractRNG,Nothing}
    )
             
    # index = max
    FIRST = true
    DECODED = false

    for i in 1:max

        if test && printing  
            println("### Iteration #$i ###")
        end

        if supermode == "FLOO"
            flooding!(d,
                      Lq,
                      Lr,
                      Lf,
                      cn2vn,
                      vn2cn,
                      Lrn,
                      signs,
                      phi
            )  
        elseif supermode == "LBP"
            LBP!(d,
                 Lr,
                 Lq,
                 Lf,
                 cn2vn,
                 vn2cn,
                 Lrn,
                 syndrome,
                 Ldn,
                 visited_vns,
                 mode == "iLBP"
            )   
        elseif supermode == "RBP"
            RBP!(
                Residues,
                d,
                Lr,
                Ms,
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
                rgn_sample
            )
            # reset factors
            Factors[H] .= 1.0
        end

        calc_syndrome!(syndrome,d,cn2vn)

        if test && printing    
                println("Max LLR estimate errors: ")
                for j in eachindex(d)
                    print(Int(d[j] != c[j]))
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
            if FIRST && iszero(syndrome)
                FIRST = false
                # index = i
                if @fastmath d == c
                    DECODED = true
                end
                if stop
                    break
                end
            end
            @fastmath bit_error .= (d .≠ c)
            @fastmath @inbounds ber[i] = sum(bit_error)
        end
    end

    # return DECODED, index
    return DECODED

end