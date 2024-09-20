################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Main loop of the SPA algorithm

include("flooding.jl")
include("llr_flooding.jl")
include("calc_syndrome.jl")
include("LBP.jl")
include("RBP.jl")

function 
    SPA!(
        d::Vector{Bool},
        ber::Vector{<:AbstractFloat},
        c::Vector{Bool},
        bit_error::Vector{Bool},
        Lr::Matrix{<:AbstractFloat}, 
        Lq::Matrix{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        nodes2checks::Vector{Vector{T}} where {T<:Integer},
        Lf::Vector{<:AbstractFloat},
        syndrome::Vector{Bool},
        Lrn::Union{Vector{<:AbstractFloat},Nothing},
        sn::Union{Vector{Int8},Nothing},        
        phi::Union{Vector{<:AbstractFloat},Nothing},
        mode::String,
        flooding::Bool,
        test::Bool,
        d_test::Union{Vector{Bool},Nothing},
        syndrome_test::Union{Vector{Bool},Nothing},
        f::Union{Matrix{<:AbstractFloat},Nothing},
        r::Union{Array{<:AbstractFloat,3},Nothing},
        q::Union{Array{<:AbstractFloat,3},Nothing},
        printing::Union{Bool,Nothing}
    )
             
    index = MAX
    FIRST = true
    DECODED = false

    for i in 1:MAX

        if flooding
            llr_flooding!(
                d,
                Lr, 
                Lq,
                checks2nodes,
                nodes2checks,
                Lf,
                Lrn,
                sn,
                phi,
            )            
        elseif mode == "LBP"
            LBP!(
                d,
                Lr,
                Lq,
                Lf,
                checks2nodes,
                nodes2checks,
                Lrn
            )
        elseif mode == "RBP"
            max_residue = 1e-16
            max_coords = [1,1]
            RBP!(
                d,
                Lr,
                max_coords,
                max_residue,
                Lq,
                Lf,
                checks2nodes,
                nodes2checks,
                Lrn
            )
        end

        calc_syndrome!(
            syndrome,
            d,
            checks2nodes
        )

        if test
            flooding!(
                d_test,
                r, 
                q,
                checks2nodes,
                nodes2checks,
                f
            )
            calc_syndrome!(
                syndrome_test,
                d_test,
                checks2nodes
            )
            if printing

                println("Iteration #$i")
    
                println("MAP SIMPLE estimate: $d_test")
    
                println("MAP Δ-LLRs estimate: $d")
    
                println("MAP SIMPLE syndrome: $syndrome_test")
    
                println("LLR Δ-LLRs syndrome: $syndrome")
    
            end
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