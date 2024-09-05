################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Simple horizontal update of the Sum-Product Algorithm

function simple_horizontal_update!(r::Array{Float64,3},
                                   δq::Matrix{Float64},
                                   indices_n::Vector{Vector{Int64}})

    M = length(indices_n)
    for m=1:M
        for n in indices_n[m]
            δr = 1
            for nn in indices_n[m]
                if nn != n
                    @inbounds δr *= δq[nn,m]
                end
            end
            @inbounds r[m,n,1] = 0.5*(1+δr)
            @inbounds r[m,n,2] = 0.5*(1-δr)
        end
    end
end