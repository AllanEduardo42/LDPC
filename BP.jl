################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Main loop of the BP routine

include("update_check2nodes_messages.jl")
include("update_node2checks_messages.jl")
include("flooding.jl")
include("calc_syndrome.jl")
include("LBP.jl")
include("iLBP.jl")
include("RBP.jl")
include("lRBP.jl")

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
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        nodes2checks::Vector{Vector{T}} where {T<:Integer},
        Lrn::Union{Vector{<:AbstractFloat},Nothing},
        sn::Union{Vector{Bool},Nothing},        
        phi::Union{Vector{<:AbstractFloat},Nothing},
        printing::Union{Bool,Nothing},        
        Residues::Union{Matrix{<:AbstractFloat},Nothing},
        Edges::Union{Matrix{<:Integer},Nothing},
        maxcoords::Union{Vector{<:Integer},Nothing},
        Factors::Union{Matrix{<:AbstractFloat},Nothing},
        pfactors::Union{AbstractFloat,Nothing},
        num_edges::Union{Integer,Nothing},
        Ldn::Union{Vector{<:AbstractFloat},Nothing},
        visited_nodes::Union{Vector{Bool},Nothing}
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
        elseif mode == "oLBP"
            LBP!(
                d,
                Lr,
                Lq,
                Lf,
                checks2nodes,
                nodes2checks,
                Lrn,
                Ldn,
                visited_nodes
            )   
        elseif mode == "iLBP"
            iLBP!(
                d,
                Lr,
                Lq,
                Lf,
                checks2nodes,
                nodes2checks,
                Lrn,
                syndrome,
                Ldn,
                visited_nodes
            )
        elseif mode == "LRBP"
            Edges .*= 0
            lRBP!(
                d,
                Lr,
                maxcoords,
                Lq,
                Lf,
                checks2nodes,
                nodes2checks,
                sn,
                Edges,
                Factors,
                pfactors,
                num_edges
            )
            # reset penalties
            reset_penalties!(Factors,checks2nodes)
        elseif mode == "oRBP"
            Edges .*= 0
            RBP!(
                d,
                Lr,
                maxcoords,
                Lq,
                Lf,
                checks2nodes,
                nodes2checks,
                sn,
                Edges,
                Factors,
                pfactors,
                num_edges,
                Residues
            )
            # reset penalties
            reset_penalties!(Factors,checks2nodes)
        end

        calc_syndrome!(
            syndrome,
            d,
            checks2nodes
        )

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

function 
    reset_penalties!(
        Factors::Matrix{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer}
    )

    check = 0
    for nodes in checks2nodes
        check += 1
        for node in nodes
            Factors[check,node] = 1.0
        end
    end
end