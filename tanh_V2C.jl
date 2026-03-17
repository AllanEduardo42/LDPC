function 
    tanh_V2C(
        arg::Float64,
        ::Float64,
        ::Float64
    )
    return arg
end

function 
    tanh_V2C(
        post_LLR::Float64,
        c2v::Float64,
        ::Nothing
    )

    # begin
    @fastmath begin
        tanh(0.5*(post_LLR - c2v))
    end
    
end
