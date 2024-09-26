################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Main loop of the BP routine

include("simple_horizontal_update.jl")
include("vertical_update_and_MAP.jl")
include("llr_horizontal_update.jl")
include("llr_vertical_update_and_MAP.jl")
include("calc_syndrome.jl")
include("LBP.jl")
include("iLBP.jl")
include("RBP.jl")
include("lRBP.jl")
include("flooding.jl")

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
        Lf::Vector{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},      
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        nodes2checks::Vector{Vector{T}} where {T<:Integer},
        Lrn::Union{Vector{<:AbstractFloat},Nothing},
        sn::Union{Vector{Bool},Nothing},        
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
        num_edges::Union{Integer,Nothing}
    )
             
    index = max
    FIRST = true
    DECODED = false

    for i in 1:max

        if mode == "FLOO"
            flooding!(
                d,
                Lq,
                Lr,
                Lf,
                checks2nodes,
                nodes2checks,
                Lrn,
                sn,
                phi
            )  
        elseif mode == "_LBP"
            LBP!(
                d,
                Lr,
                Lq,
                Lf,
                checks2nodes,
                nodes2checks,
                Lrn
            )   
        elseif mode == "ILBP"
            iLBP!(
                d,
                Lr,
                Lq,
                Lf,
                checks2nodes,
                nodes2checks,
                Lrn,
                syndrome
            )
        elseif mode == "LRBP"
            Edges .*= 0
            lRBP!(
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
                penalty_factor,
                num_edges
            )
            # reset penalties
            reset_penalties!(penalty,checks2nodes)
        elseif mode == "_RBP"
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
                penalty_factor,
                num_edges,
                R
            )
            # reset penalties
            reset_penalties!(penalty,checks2nodes)
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

function 
    reset_penalties!(
        penalty::Matrix{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer}
    )

    check = 0
    for nodes in checks2nodes
        check += 1
        for node in nodes
            penalty[check,node] = 1.0
        end
    end
end