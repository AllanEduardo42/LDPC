################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Vertical update and MAP estimate of the LLR based Sum-Product Algorithm (with
# "Inf" restriction)

function
    update_node2checks_messages!(
        q::Array{<:AbstractFloat,3},
        r::Array{<:AbstractFloat,3},
        Ld::Vector{<:AbstractFloat},
        node::Integer,
        checks::Vector{<:Integer}
    )

    
    for check in checks
        @inbounds @fastmath Ld[1] *= r[check,node,1]
        @inbounds @fastmath Ld[2] *= r[check,node,2]
    end
    for check in checks
        @inbounds @fastmath q1 = Ld[1] / (r[check,node,1] + eps())  
        @inbounds @fastmath q2 = Ld[2] / (r[check,node,2] + eps())  
        @fastmath a = q1 + q2
        @inbounds @fastmath q[check,node,1] = q1/a
        @inbounds @fastmath q[check,node,2] = q2/a
    end

    return signbit(Ld[1]-Ld[2])

end

function
    update_node2checks_messages!(
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lf::AbstractFloat,
        node::Integer,
        checks::Vector{<:Integer}
    )

    Ld = calc_Ld(
        node,
        checks,
        Lf,
        Lr
    )
    for check in checks
        @inbounds @fastmath Lq[check,node] = Ld - Lr[check,node]
    end

    return Ld, signbit(Ld)
end

function
    update_node2check_message(
        nodes2checks_nmax::Vector{<:Integer},
        node::Integer,
        check::Integer,
        Lq::AbstractFloat,
        Lr::Matrix{<:AbstractFloat}
    )
    
    for c in nodes2checks_nmax
        if c != check
            @inbounds @fastmath Lq += Lr[c,node]
        end
    end

    return Lq
end

function 
    calc_Ld(
        node::Integer,
        checks::Vector{<:Integer},
        Lf::AbstractFloat,
        Lr::Matrix{<:AbstractFloat}
    )
    Ld = Lf
    for check in checks
        @inbounds @fastmath Ld += Lr[check,node]
    end
    
    return Ld

end

function 
    MAP!(
        d::Vector{Bool},
        nodes2checks::Vector{Vector{T}} where {T<:Integer},
        Lf::Vector{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat}
    )
    
    node = 0
    Ld = 0.0
    for checks in nodes2checks
        node += 1
        Ld = calc_Ld(
                    node,
                    checks,
                    Lf[node],
                    Lr
                )
        @inbounds d[node] = signbit(Ld)
    end

end

