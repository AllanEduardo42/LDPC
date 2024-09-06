################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# MAin loop of the SPA algorithm

function SPA!(d::Vector{Int64},
              ber::Vector{Float64},
              c::Vector{Int64},
              x::BitArray,
              Lr::Matrix{Float64}, 
              Lq::Matrix{Float64},
              indices_n::Vector{Vector{Int64}},
              indices_m::Vector{Vector{Int64}},
              ΔLf::Vector{Float64},
              syndrome::Vector{Int64},
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

        calc_syndrome!(syndrome,indices_n,d)
        
        s = sum(syndrome)

        if FIRST && s == 0
            FIRST = false
            index = i
            if d == c
                DECODED = true
            end
        end
        x .= d .!= c
        ber[i] = sum(x)

    end

    return DECODED, index

end

################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# MAin loop of the SPA algorithm

function SPA!(d::Vector{Int64},
              ber::Vector{Float64},
              c::Vector{Int64},
              x::BitArray,
              Lr::SparseMatrixCSC{Float64, Int64}, 
              Lq::SparseMatrixCSC{Float64, Int64},
              indices_n::Vector{Vector{Int64}},
              indices_m::Vector{Vector{Int64}},
              ΔLf::Vector{Float64},
              syndrome::Vector{Int64},
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

        calc_syndrome!(syndrome,indices_n,d)
        
        s = sum(syndrome)

        if FIRST && s == 0
            FIRST = false
            index = i
            if d == c
                DECODED = true
            end
        end
        x .= d .!= c
        ber[i] = sum(x)

    end

    return DECODED, index

end

