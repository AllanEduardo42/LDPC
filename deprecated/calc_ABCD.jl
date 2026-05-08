################ ALTERNATIVE TO HYPERBOLIC TANGENT USING TABLE ################
function calc_ABCD!(
    V2C::Matrix{Float64},
    ci::Int,
    Nci::Vector{Int},
    signs::Vector{Bool},
    phi::Vector{Float64}
)

    sum_c2v = 0.0
    s = false
    @fastmath @inbounds for vj in Nci
        v2c = V2C[ci,vj]
        sig = signbit(v2c)
        s ⊻= sig
        signs[vj] = sig
        sum_c2v += ϕ(abs(v2c),phi) 
    end

    return sum_c2v, s, nothing, nothing
end

function calc_C2V(
    sum_c2v::Float64,   #A
    s::Bool,            #B          
    ::Nothing,          #C
    ::Nothing,          #D
    vj::Int,
    v2c::Float64,
    signs::Vector{Bool},
    phi::Vector{Float64}
)

    @fastmath @inbounds begin
        x = abs(sum_c2v - v2c)
        y = signs[vj] ⊻ s
        return (1 - 2*y)*ϕ(x,phi)
    end
end