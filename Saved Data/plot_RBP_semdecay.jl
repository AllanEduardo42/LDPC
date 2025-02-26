using DelimitedFiles
using Plots

SNR = [1.2, 1.6, 1.8, 2.0]

p1 = readdlm("./Saved Data/RBP absoluto sem decay/FER_RBP 1.0.txt",'\t',Float64,'\n')
p2 = readdlm("./Saved Data/RBP relativo sem decay/FER_RBP 1.0.txt",'\t',Float64,'\n')

plotlyjs()
lim = log10(1/maximum(1000*2^10))
labels = Vector{String}()
for snr in SNR
    push!(labels,"SNR (dB) = $snr, absolute")
end
labels = permutedims(labels)
titlefer = "FER (decay factor 1.0)"
plot(
    1:20,
    p1,
    xlabel="Iteration",
    label=labels,
    lw=2,
    title=titlefer,
    ylims=(lim,0)
)
labels = Vector{String}()
for snr in SNR
    push!(labels,"SNR (dB) = $snr, relative")
end
labels = permutedims(labels)
titlefer = "FER RBP (decay factor 1.0)"
plot!(
    1:20,
    p2,
    xlabel="Iteration",
    label=labels,
    lw=2,
    ls=:dot,
    title=titlefer,
    ylims=(lim,0),
    color=[1 2 3 4]
)