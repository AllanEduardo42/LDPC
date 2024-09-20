################################################################################
# Allan Eduardo Feitosa
# 19 set 2024
# Core of the horizontal update of the LLR based Sum-Product Algorithm 
# using the tang function (with Inf restriction)
#(Obs: this core is a separate function because it is also used in the LBP
# algorithm)

function
    tanh_llr_horizontal_update!(
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        nodes::Vector{<:Integer},
        check::Integer,
        Lrn::Vector{<:AbstractFloat},
    )
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