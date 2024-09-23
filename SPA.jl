################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Main loop of the SPA algorithm

include("simple_horizontal_update.jl")
include("vertical_update_and_MAP.jl")
include("llr_horizontal_update.jl")
include("llr_vertical_update_and_MAP.jl")
include("calc_syndrome.jl")
include("LBP.jl")
include("RBP.jl")
include("RBP_R.jl")

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
        printing::Union{Bool,Nothing},
        max::Integer,
        R::Union{Matrix{<:AbstractFloat},Nothing}
    )
             
    index = max
    FIRST = true
    DECODED = false

    if mode == "RBP_R"
        R .*= 0.0
    end

    for i in 1:max

        if flooding
            llr_horizontal_update!(
                Lr,
                Lq,
                checks2nodes,
                Lrn,
                sn,
                phi
            )
            llr_vertical_update_and_MAP!(
                Lq,
                d,
                Lr,
                Lf,
                nodes2checks
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
            max_coords = [1,1]
            RBP!(
                d,
                Lr,
                max_coords,
                Lq,
                Lf,
                checks2nodes,
                nodes2checks,
                sn
            )
        elseif mode == "RBP_R"
            max_coords = [1,1]
            RBP_R!(
                d,
                Lr,
                max_coords,
                Lq,
                Lf,
                checks2nodes,
                nodes2checks,
                sn,
                R
            )
        end

        calc_syndrome!(
            syndrome,
            d,
            checks2nodes
        )

        if test            
            ### Conventional simplified SPA
            δQ = q[:,:,1]-q[:,:,2]    
            simple_horizontal_update!(
                r,
                δQ,
                checks2nodes
            )
            vertical_update_and_MAP!(
                q,
                d_test,
                r,
                f,
                nodes2checks
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