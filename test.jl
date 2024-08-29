function test(x)
    s = 0
    for i in x
        s = i + 1
    end
end

function test2(x)
    s = 0
    for i in eachindex(x)
        s = i + x[i]
    end
end