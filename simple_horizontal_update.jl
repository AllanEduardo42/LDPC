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

    check = 0
    for nodes in checks2nodes
        check += 1
        for node in nodes
            δr = 1
            for n in nodes
                if n != node
                    @inbounds δr *= δq[check,n]
                end
            end
            @inbounds r[check,node,1] = 0.5*(1+δr)
            @inbounds r[check,node,2] = 0.5*(1-δr)
        end
    end
end