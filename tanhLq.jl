function 
    tanhLq(
        Lq::Float64,
        ::Vector{Bool}
    )
    return Lq
end

function 
    tanhLq(
        Lq::Float64,
        ::Nothing
    )
    return tanh(0.5*Lq)
end
