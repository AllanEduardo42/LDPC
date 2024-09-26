function f2(x,s,arg)
    y = zeros(length(x))
    for i in eachindex(x)
        y[i] = (2*(x[i] == s)-1)*arg
    end
end

function f1(x,s,arg)
    y = zeros(length(x))
    for i in eachindex(x)
        y[i] = (1 - 2*(x[i] ⊻ s))*arg
    end
end

N = 2000
t1 = zeros(N)
t2 = zeros(N)

for i in 1:N
    x = @timed f1(rand(Bool,1000000),false,3.14)
    t1[i] = x.time
    x = @timed f2(rand(Bool,1000000),false,3.14)
    t2[i] = x.time
end

histogram(
    [t1 t2],
    layout=grid(2,1),
    xlim=(0.001,0.002)
)