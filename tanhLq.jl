function 
    tanhLq(
        arg::Float64,
        ::Float64,
        ::Float64
    )
    return arg
end

function 
    tanhLq(
        ld::Float64,
        lr::Float64,
        ::Nothing
    )

    # begin
    @fastmath begin
        tanh(0.5*(ld - lr))
    end
    
end
