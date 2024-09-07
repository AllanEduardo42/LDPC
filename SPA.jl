################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# MAin loop of the SPA algorithm

function SPA!(d::Vector{Bool},
              ber::Vector{Float64},
              c::Vector{Bool},
              bit_error::Vector{Bool},
              Lr::Matrix{Float64}, 
              Lq::Matrix{Float64},
              indices_n::Vector{Vector{Int64}},
              indices_m::Vector{Vector{Int64}},
              ΔLf::Vector{Float64},
              syndrome::Vector{Bool},
              sn::Any,
              Lrn::Any,
              phi::Any)
    
    # varargs = (sn::Vector{Int64},
    #            Lrn::Vector{Float64},
    #            phi::phi::Vector{Float64})
             
    index = MAX
    FIRST = true
    DECODED = false

    for i in 1:MAX
        
        llr_horizontal_update!(Lr,Lq,indices_n,sn,Lrn,phi)
        llr_vertical_update_and_MAP!(Lq,d,Lr,ΔLf,indices_m)

        calc_syndrome!(syndrome,d,indices_n)

        if FIRST && iszero(syndrome)
            FIRST = false
            index = i
            if d == c
                DECODED = true
            end
        end
        bit_error .= (d .!= c)
        ber[i] = sum(bit_error)

    end

    return DECODED, index

end