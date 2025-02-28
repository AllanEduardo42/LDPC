using DelimitedFiles
using Plots

SNR = [1.2, 1.6, 1.8, 2.0]
plotlyjs()
lim = log10(1/maximum(1000*2^10))

p = readdlm("./Saved Data/RBP absoluto e relativo/FER_RBP 0.8 absoluto.txt",'\t',Float64,'\n')
labels = Vector{String}()
for snr in SNR
    push!(labels,"SNR (dB) = $snr, absolute")
end
labels = permutedims(labels)
titlefer = "FER RBP (decay factor 0.8)"
plot(
    1:20,
    p,
    xlabel="Iteration",
    label=labels,
    lw=2,
    title=titlefer,
    ylims=(lim,0)
)

p = readdlm("./Saved Data/RBP absoluto e relativo/FER_RBP 0.8 relativo.txt",'\t',Float64,'\n')
labels = Vector{String}()
for snr in SNR
    push!(labels,"SNR (dB) = $snr, relative")
end
labels = permutedims(labels)
titlefer = "FER (decay factor 0.8)"
display(plot!(
    1:20,
    p,
    xlabel="Iteration",
    label=labels,
    lw=2,
    ls=:dot,
    title=titlefer,
    ylims=(lim,0),
    color=[1 2 3 4]
))

p = readdlm("./Saved Data/RBP absoluto e relativo/FER_RBP 1.0 absoluto.txt",'\t',Float64,'\n')
labels = Vector{String}()
for snr in SNR
    push!(labels,"SNR (dB) = $snr, absolute")
end
labels = permutedims(labels)
titlefer = "FER RBP (decay factor 1.0)"
plot(
    1:20,
    p,
    xlabel="Iteration",
    label=labels,
    lw=2,
    title=titlefer,
    ylims=(lim,0)
)

p = readdlm("./Saved Data/RBP absoluto e relativo/FER_RBP 1.0 relativo.txt",'\t',Float64,'\n')
labels = Vector{String}()
for snr in SNR
    push!(labels,"SNR (dB) = $snr, relative")
end
labels = permutedims(labels)
titlefer = "FER (decay factor 1.0)"
plot!(
    1:20,
    p,
    xlabel="Iteration",
    label=labels,
    lw=2,
    ls=:dot,
    title=titlefer,
    ylims=(lim,0),
    color=[1 2 3 4]
)