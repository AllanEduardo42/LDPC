x = zeros(1000000)
y = zeros(1000000)

for i in 1:1000000
    x[i] = @allocated trues(i)
    y[i] = @allocated ones(Bool,i)
end

plot(x, label="trues")
plot!(y, labels="bools")