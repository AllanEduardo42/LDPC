using LinearAlgebra
using Random
using Statistics
using Plots
using Distributions
using StatsPlots

N = 10000
v = 2
m = -1

z = zeros(N)

u = zeros(N)

max_checks = 6

plotlyjs()

p = plot(layout=(6,max_checks),size= 4 .*(600,500))

for vars in 1:6

    x = zeros(N,vars)

    y = zeros(N,vars)

    global u *= 0

    for checks = 1:max_checks

        x = m .+ sqrt(v)*randn(N,vars)

        xx = 2*x/v

        y = tanh.(xx/2)

        global z = prod(y,dims=2)

        global u += 2*atanh.(z)

        d = Distributions.fit(Normal,u)

        mea = d.μ

        sig = d.σ

        max = 1/sqrt(2*pi*sig^2)

        histogram!(
            u,
            normalize=true,
            title="var=$vars, chk= $(checks)",
            subplot = (vars-1)*max_checks + checks,
            legend=false)

        plot!([mea, mea], [0, max], subplot = (vars-1)*max_checks + checks,lw=3)
        plot!(d,subplot = (vars-1)*max_checks + checks,lw=3)
        

    end

end

display(p)

# plot(Tuple(P),layout=(2,3))

# display(p)

# U = mean(u)

## min-sum

# s = sign.(x)

# c = prod(s,dims=2)

# v = c .* minimum(abs.(x), dims=2)

# t = -1:0.01:3

# f = 1/sqrt(2*pi*v) * exp.( -(t .- m).^2/(2*v))

# p = plot!(t,f,lw=2)

# p = plot!([U, U],[0,8000], lw=2)

# display(p)

# p = histogram!(v)

# display(p)

# e = 1 ./(1 .+ exp.(-u))

# p = histogram(e)

# display(p)

# p = histogram(z)

# display(p)


