using Plots

N = 100000

x = rand([0,1],N)
σ = 0.8
n = σ*randn(N)
tx = 2*x .- 1 .+ n

f = zeros(N,2)
k = 1/(sqrt(2π)*σ)
s = 2σ^2
for i in eachindex(tx)
    f[i,1] = k*exp(-(tx[i]+1)^2/s)
    f[i,2] = k*exp(-(tx[i]-1)^2/s)
end

L1 = -log.(f[:,1])
L2 = -log.(f[:,2])
ΔL = abs.(L1 - L2)

f_minus = abs.(log.(1 .- exp.(-L1)))
f_plus = abs.(log.(1 .+ exp.(-L1)))
Δf_minus = abs.(log.(1 .- exp.(-ΔL)))
Δf_plus = abs.(log.(1 .+ exp.(-ΔL)))

plotlyjs()

display(histogram(f_minus))
display(histogram(f_plus))
display(histogram(Δf_minus))
histogram(Δf_plus)

