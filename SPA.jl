################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# MAin loop of the SPA algorithm

function SPA(S::Int64,
             i::Int64,
             M::Int64,
             N::Int64,
             Lr::Matrix{Float64}, 
             Lq::Matrix{Float64},
             indices_n::Vector{Vector{Int64}},
             indices_m::Vector{Vector{Int64}},
             d::Vector{Int64},
             ΔLf::Vector{Float64},
             syndrome::Vector{Int64},
             varargs...)
    
    # varargs = (sn::Vector{Int64},
    #            Lrn::Vector{Float64},
    #            phi::phi::Vector{Float64})
             

    while S != 0 && i < MAX

        i += 1

        Lr = llr_horizontal_update(M,Lr,Lq,indices_n,varargs...)
        Lq, d = llr_vertical_update_and_MAP(Lq, Lr, d, N, ΔLf, indices_m)

        syndrome = calc_syndrome(syndrome,M,indices_n,d)
        
        S = sum(syndrome)

    end

    return S, i

end

function calc_syndrome(syndrome::Vector{Int64},
                       M::Int64,
                       indices_n::Vector{Vector{Int64}},
                       d::Vector{Int64})

    syndrome .*= 0
    for m=1:M
        for n in indices_n[m]
            syndrome[m] += d[n]
        end
    end

    syndrome .%= 2 

    return syndrome

end
