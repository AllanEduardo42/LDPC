################################################################################
# Allan Eduardo Feitosa
# 17 set 2024
# LBP Sum-Product Algorithm

function
    LBP!(
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        d::Vector{Bool},
        ΔLf::Vector{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        nodes2checks::Vector{Vector{T}} where {T<:Integer},
        Lrn::Vector{<:AbstractFloat},
    )

    check = 0
    for nodes in checks2nodes
        check += 1
        # vertical update        
        for node in nodes
            @inbounds Lq[check,node] = ΔLf[node]
            for c in nodes2checks[node]
                if c != check
                    @inbounds @fastmath Lq[check,node] += Lr[c,node]
                end
            end
        end
        # horizontal update
        pLr = 1.0
        for node in nodes
            @inbounds @fastmath Lrn[node] = tanh(0.5*Lq[check,node])
            @inbounds @fastmath pLr *= Lrn[node]
        end
        for node in nodes
            @inbounds @fastmath x = pLr/Lrn[node]
            if abs(x) < 1
                @inbounds @fastmath Lr[check,node] = 2*atanh(x)
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
    SPA_LBP!(
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