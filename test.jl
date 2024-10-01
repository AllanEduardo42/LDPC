using Plots

function f1!(x,y)    
    for i in eachindex(x)
        y[i] = log(1 + 2/(exp(i) - 1))
    end
end

function f2!(x,y)
    for i in eachindex(x)
        m = exp(i) - 1
        y[i] = log(1 + 2/m)
    end
end

N = 1000
L = 1000000
t1 = zeros(N)
t2 = zeros(N)
y = Vector{Float64}(undef,L)

for i in 1:N
    x = randn(L)
    t = @timed f1!(x,y)
    t1[i] = t.time
end

for i in 1:N
    x = randn(L)
    t = @timed f2!(x,y)
    t2[i] = t.time
end

plotlyjs()

histogram(
    [t1 t2],
    layout=grid(2,1),
    xlim=(0.01,0.011),
    ylim=(0, 600)
)