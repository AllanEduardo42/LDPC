################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Main loop of the BP routine

include("auxiliary_functions.jl")
include("update_Lr.jl")
include("update_Lq.jl")
include("flooding.jl")
include("calc_syndrome.jl")
include("LBP.jl")
include("iLBP.jl")
include("RBP.jl")
include("LRBP.jl")

function 
    BP!(
        mode::String,
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
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        Lrn::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},        
        phi::Union{Vector{<:AbstractFloat},Nothing},
        printing::Union{Bool,Nothing},        
        Residues::Union{Matrix{<:AbstractFloat},Nothing},
        Edges::Union{Matrix{<:Integer},Nothing},
        maxcoords::Union{Vector{<:Integer},Nothing},
        Factors::Union{Matrix{<:AbstractFloat},Nothing},
        pfactors::Union{AbstractFloat,Nothing},
        num_edges::Union{Integer,Nothing},
        Ldn::Union{Vector{<:AbstractFloat},Nothing},
        visited_vns::Union{Vector{Bool},Nothing},
        samples::Union{Vector{<:Integer},Nothing}
    )
             
    index = max
    FIRST = true
    DECODED = false

    for i in 1:max

        if mode == "FLOO"
            flooding!(d,Lq,Lr,Lf,cn2vn,vn2cn,Lrn,signs,phi)  
        elseif mode == "LBP"
            LBP!(d,Lr,Lq,Lf,cn2vn,vn2cn,Lrn,Ldn,visited_vns)   
        elseif mode == "iLBP"
            iLBP!(d,Lr,Lq,Lf,cn2vn,vn2cn,Lrn,syndrome,Ldn,visited_vns)
        elseif mode == "RBP"
            Edges .*= 0
            RBP!(
                d,
                Lr,
                maxcoords,
                Lq,
                Lf,
                cn2vn,
                vn2cn,
                signs,
                Edges,
                Factors,
                pfactors,
                num_edges,
                Residues,
                samples
            )
            # reset factors
            reset_factors!(Factors,cn2vn)
        elseif mode == "LRBP"
            Edges .*= 0
            LRBP!(
                d,
                Lr,
                maxcoords,
                Lq,
                Lf,
                cn2vn,
                vn2cn,
                signs,
                Edges,
                Factors,
                pfactors,
                num_edges,
                Ldn,
                syndrome
            )
            # reset factors
            reset_factors!(Factors,cn2vn)
        end

        calc_syndrome!(syndrome,d,cn2vn)

        if test && printing
                println("Iteration #$i")     
                println("MAP estimate: $d")      
                println("Syndrome: $syndrome")    
        else
            if FIRST && iszero(syndrome)
                FIRST = false
                index = i
                if d == c
                    DECODED = true
                end
            end
            bit_error .= (d .≠ c)
            ber[i] = sum(bit_error)
        end
    end

    return DECODED, index

end