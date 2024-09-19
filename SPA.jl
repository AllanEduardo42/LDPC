################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# MAin loop of the SPA algorithm

function 
    SPA!(
        d::Vector{Bool},
        ber::Vector{<:AbstractFloat},
        c::Vector{Bool},
        bit_error::Vector{Bool},
        Lr::Matrix{<:AbstractFloat}, 
        Lq::Matrix{<:AbstractFloat},
        indices_row::Vector{Vector{T}} where {T<:Integer},
        indices_col::Vector{Vector{T}} where {T<:Integer},
        ΔLf::Vector{<:AbstractFloat},
        syndrome::Vector{Bool},
        sn::Union{Vector{Int8},Nothing},
        Lrn::Union{Vector{<:AbstractFloat},Nothing},
        phi::Union{Vector{<:AbstractFloat},Nothing}
    )
    
    # varargs = (sn::Vector{<:Integer},
    #            Lrn::Vector{<:AbstractFloat},
    #            phi::phi::Vector{<:AbstractFloat})
             
    index = MAX
    FIRST = true
    DECODED = false

    for i in 1:MAX
        
        llr_horizontal_update!(
            Lr,
            Lq,
            indices_row,
            sn,
            Lrn,
            phi
        )
        llr_vertical_update_and_MAP!(
            Lq,
            d,
            Lr,
            ΔLf,
            indices_col
        )

        calc_syndrome!(
            syndrome,
            d,
            indices_row
        )

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

    return DECODED, index

end