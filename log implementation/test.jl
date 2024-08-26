using Plots

N = 10000

x = rand(N)
y = rand(N)
Lx = -log.(x)
Ly = -log.(y)
ΔL = abs.(Lx - Ly)

range = 1/N:1/N:10
f = exp.(-range)

plotlyjs()

histogram(Lx, normalize=:pdf)
display(plot!(range, f, linewidth=2))
histogram(ΔL, normalize=:pdf)
plot!(range, f, linewidth=2)

