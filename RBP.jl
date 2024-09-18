################################################################################
# Allan Eduardo Feitosa
# 17 set 2024
# RBP Sum-Product Algorithm

function
    RBP!(
        Lr::Matrix{Float64},
        max_coords::Vector{Int64},
        max_residue::Float64,
        Lq::Matrix{Float64},
        d::Vector{Bool},
        ΔLf::Vector{Float64},
        checks2nodes::Vector{Vector{Int64}},
        nodes2checks::Vector{Vector{Int64}},
        Lrn::Vector{Float64},
    )

    for m in eachindex(checks2nodes)

        cmax = max_coords[1]
        nmax = max_coords[2]

        # update the message with largest residue
        @inbounds @fastmath Lr[cmax,nmax] += max_residue
        max_residue = 0
        
        for check in nodes2checks[nmax]
            if check ≠ cmax
                # vertical update of Lq[check,nmax], check ≠ cmax
                @inbounds Lq[check,nmax] = ΔLf[nmax]
                for c in nodes2checks[nmax]
                    if c != check
                        @inbounds @fastmath Lq[check,nmax] += Lr[c,nmax]
                    end
                end
                # horizontal update of Lr[node,check], node ≠ nmax, check ≠ cmax
                pLr = 1.0
                for node in checks2nodes[check]
                    @inbounds @fastmath Lrn[node] = tanh(0.5*Lq[check,node])
                    @inbounds @fastmath pLr *= Lrn[node]
                end
                for node in checks2nodes[check]
                    if node != nmax
                        @inbounds @fastmath L = 2*atanh(pLr/Lrn[node])
                        @inbounds @fastmath ΔLr = L - Lr[check,node]
                        if abs(ΔLr) > abs(max_residue)
                            max_residue = ΔLr
                            max_coords[1] = check
                            max_coords[2] = node
                        end
                    end
                end
            end
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
    SPA_RBP!(
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

    max_coords = ones(Int,2)
    max_residue = 0.0

    for i in 1:MAX
        
        RBP!(
            Lr,
            max_coords,
            max_residue,
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