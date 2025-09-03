function 
    tanhLq(
        arg::Float64,
        ::Float64
    )
    return arg
end

function 
    tanhLq(
        arg::Float64,
        ::Nothing
    )
    return tanh(0.5*arg)
end
