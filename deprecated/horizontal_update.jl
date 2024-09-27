################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Horizontal update of the Sum-Product Algorithm

function
    horizontal_update!(
        r::Array{<:AbstractFloat,3},
        q::Array{<:AbstractFloat,3},
        indices_row::Vector{Vector{T}} where {T<:Integer}
    )

    m = 0
    for indices in indices_row
        m += 1
        S = length(indices)-1
        for n in indices           
            for s = 0:2^S-1
                dig = digits(s, base = 2, pad = S)
                count = 0
                rr = 1
                if iseven(sum(dig))
                    for nn in indices
                        if nn != n
                            count += 1
                            @inbounds @fastmath rr *= q[nn,m,dig[count]+1]
                        end
                    end
                    @inbounds @fastmath r[m,n,1] += rr
                else
                    for nn in indices
                        if nn != n
                            count += 1
                            @inbounds @fastmath rr *= q[nn,m,dig[count]+1]
                        end               
                    end
                    @inbounds @fastmath r[m,n,2] += rr
                end
            end
        end
    end

end