function f1(x)
    y = zeros(length(x))
    for i in eachindex(x)
        y[i] = log(1 + 2/(exp(i) - 1))
    end
end

function f2(x)
    y = zeros(length(x))
    for i in eachindex(x)
        m = exp(i) - 1
        y[i] = log(1 + 2/m)
    end
end

N = 1000
t1 = zeros(N)
t2 = zeros(N)

for i in 1:N
    x = @timed f1(rand(1000000))
    t1[i] = x.time
end

for i in 1:N
    x = @timed f2(rand(1000000))
    t2[i] = x.time
end

histogram(
    [t1 t2],
    layout=grid(2,1),
    xlim=(0.01,0.015),
    ylim=(0, 600)
)