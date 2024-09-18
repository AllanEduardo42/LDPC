################################################################################
# Allan Eduardo Feitosa
# 17 set 2024
# LBP Sum-Product Algorithm

function
    LBP!(
        Lr::Matrix{Float64},
        Lq::Matrix{Float64},
        d::Vector{Bool},
        ΔLf::Vector{Float64},
        checks2nodes::Vector{Vector{Int64}},
        nodes2checks::Vector{Vector{Int64}},
        Lrn::Vector{Float64},
    )

    for check in eachindex(checks2nodes)
        # vertical update        
        for node in checks2nodes[check]
            @inbounds Lq[check,node] = ΔLf[node]
            for c in nodes2checks[node]
                if c != check
                    @inbounds @fastmath Lq[check,node] += Lr[c,node]
                end
            end
        end
        # horizontal update
        pLr = 1.0
        for node in checks2nodes[check]
            @inbounds @fastmath Lrn[node] = tanh(0.5*Lq[check,node])
            @inbounds @fastmath pLr *= Lrn[node]
        end
        for node in checks2nodes[check]
            @inbounds @fastmath Lr[check,node] = 2*atanh(pLr/Lrn[node])
        end
    end

    #MAP
    d .*= false
    for node in eachindex(nodes2checks)
        @inbounds Ld = ΔLf[node]
        for check in nodes2checks[node]
            @inbounds @fastmath Ld += Lr[check,node]
        end
        if Ld < 0
            @inbounds d[node] = 1
        end
    end
end

function 
    SPA_LBP!(
        d::Vector{Bool},
        ber::Vector{Float64},
        c::Vector{Bool},
        bit_error::Vector{Bool},
        Lr::Matrix{Float64}, 
        Lq::Matrix{Float64},
        checks2nodes::Vector{Vector{Int64}},
        nodes2checks::Vector{Vector{Int64}},
        ΔLf::Vector{Float64},
        syndrome::Vector{Bool},
        Lrn::Vector{Float64}
    )
    
    # varargs = (sn::Vector{Int64},
    #            Lrn::Vector{Float64},
    #            phi::phi::Vector{Float64})
             
    index = MAX
    FIRST = true
    DECODED = false

    for i in 1:MAX
        
        LBP!(
            Lr,
            Lq,
            d,
            ΔLf,
            checks2nodes,
            nodes2checks,
            Lrn
        )

        calc_syndrome!(
            syndrome,
            d,
            checks2nodes
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