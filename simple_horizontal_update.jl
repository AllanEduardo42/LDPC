################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Simple horizontal update of the Sum-Product Algorithm

function 
    simple_horizontal_update!(
        r::Array{<:AbstractFloat,3},
        δq::Matrix{<:AbstractFloat},
        checks2nodes::Vector{Vector{T}} where {T<:Integer}
    )

    m = 0
    for indices in checks2nodes
        m += 1
        for n in indices
            δr = 1
            for nn in indices
                if nn != n
                    @inbounds δr *= δq[m,nn]
                end
            end
            @inbounds r[m,n,1] = 0.5*(1+δr)
            @inbounds r[m,n,2] = 0.5*(1-δr)
        end
    end
end