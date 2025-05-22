E = sum(HH)

degrees_vn = sum(HH,dims=1)
degrees_cn = sum(HH,dims=2)

λ = zeros(maximum(degrees_vn))
ρ = zeros(maximum(degrees_cn))

for i in eachindex(λ)
    λ[i] = sum(degrees_vn .== i)*i/E
end

for i in eachindex(ρ)
    ρ[i] = sum(degrees_cn .== i)*i/E
end

display([λ[end:-1:1] LAMBDA])
display(ρ[end:-1:1])