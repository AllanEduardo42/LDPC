################################################################################
# Allan Eduardo Feitosa
# 17 set 2024
# RBP Sum-Product Algorithm

function
    RBP!(
        Lr::Matrix{<:AbstractFloat},
        max_coords::Vector{<:Integer},
        max_residue::AbstractFloat,
        Lq::Matrix{<:AbstractFloat},
        d::Vector{Bool},
        ΔLf::Vector{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        nodes2checks::Vector{Vector{T}} where {T<:Integer},
        Lrn::Vector{<:AbstractFloat},
    )

    for m in eachindex(checks2nodes)

        if max_residue == 0
            break
        end

        (cmax,nmax) = max_coords

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
                # find max ΔLr[node,check], node ≠ nmax, check ≠ cmax
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
        ber::Vector{<:AbstractFloat},
        c::Vector{Bool},
        bit_error::Vector{Bool},
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        nodes2checks::Vector{Vector{T}} where {T<:Integer},
        ΔLf::Vector{<:AbstractFloat},
        syndrome::Vector{Bool},
        Lrn::Vector{<:AbstractFloat}
    )
    
    # varargs = (sn::Vector{<:Integer},
    #            Lrn::Vector{<:AbstractFloat},
    #            phi::phi::Vector{<:AbstractFloat})
             
    index = MAX
    FIRST = true
    DECODED = false

    for i in 1:MAX

        max_residue = 1e-16
        max_coords = [1,1]
        
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