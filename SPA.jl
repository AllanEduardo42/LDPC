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
        mode::String,
        test::Bool,
        max::Integer,
        syndrome::Vector{Bool},
        d::Vector{Bool},
        c::Vector{Bool},
        bit_error::Vector{Bool},
        ber::Vector{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},      
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        nodes2checks::Vector{Vector{T}} where {T<:Integer},
        Lrn::Union{Vector{<:AbstractFloat},Nothing},
        sn::Union{Vector{Int8},Nothing},        
        phi::Union{Vector{<:AbstractFloat},Nothing},
        f::Union{Matrix{<:AbstractFloat},Nothing},
        q::Union{Array{<:AbstractFloat,3},Nothing},
        r::Union{Array{<:AbstractFloat,3},Nothing},
        printing::Union{Bool,Nothing},        
        R::Union{Matrix{<:AbstractFloat},Nothing},
        Edges::Union{Matrix{<:Integer},Nothing},
        max_coords::Union{Vector{<:Integer},Nothing},
        penalty::Union{Matrix{<:AbstractFloat},Nothing},
        penalty_factor::Union{AbstractFloat,Nothing},
    )
             
    index = max
    FIRST = true
    DECODED = false

    for i in 1:max

        if mode == "TNH"
            llr_horizontal_update_tnh!(
                Lr,
                Lq,
                checks2nodes,
                Lrn
            )
            llr_vertical_update_and_MAP!(
                Lq,
                d,
                Lr,
                Lf,
                nodes2checks
            )  
        elseif mode == "ALT"
            llr_horizontal_update_alt!(
                Lr,
                Lq,
                checks2nodes,
                Lrn,
                sn
            )
            llr_vertical_update_and_MAP!(
                Lq,
                d,
                Lr,
                Lf,
                nodes2checks
            )  
        elseif mode == "TAB"
            llr_horizontal_update_tab!(
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
        elseif mode == "MIN"
            min_sum!(
                Lr,
                Lq,
                checks2nodes,
                sn
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
            Edges .*= 0
            RBP!(
                d,
                Lr,
                max_coords,
                Lq,
                Lf,
                checks2nodes,
                nodes2checks,
                sn,
                Edges,
                penalty,
                penalty_factor
            )
        elseif mode == "RBP_R"
            Edges .*= 0
            RBP_R!(
                d,
                Lr,
                max_coords,
                Lq,
                Lf,
                checks2nodes,
                nodes2checks,
                sn,
                R,
                Edges
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
            d_test = vertical_update_and_MAP!(
                q,
                r,
                f,
                nodes2checks
            )
            syndrome_test = calc_syndrome(
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