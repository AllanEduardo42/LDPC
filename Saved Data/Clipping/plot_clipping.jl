using DelimitedFiles
using Plots

SNR = collect(1.2:0.4:2.0)
plotlyjs()
lim = log10(1/maximum(1000*2^10))

p = readdlm("./Saved Data/Clipping/FER_List-RBP 0.8.txt",'\t',Float64,'\n')
labels = Vector{String}()
for snr in SNR
    push!(labels,"SNR = $snr, listsizes = 16,1")
end
labels = permutedims(labels)
titlefer = "FER (decay factor 0.8)"
plot(
    1:20,
    p,
    xlabel="Iteration",
    label=labels,
    lw=2,
    title=titlefer,
    ylims=(lim,0)
)
p = readdlm("./Saved Data/Clipping/FER_RBP 0.8.txt",'\t',Float64,'\n')
labels = Vector{String}()
for snr in SNR
    push!(labels,"SNR = $snr, no list")
end
labels = permutedims(labels)
plot!(
    1:20,
    p,
    xlabel="Iteration",
    label=labels,
    lw=2,
    ls=:dot,
    title=titlefer,
    ylims=(lim,0)
)

p = readdlm("./Saved Data/Clipping/FER_Local-RBP 0.8.txt",'\t',Float64,'\n')
labels = Vector{String}()
for snr in SNR
    push!(labels,"SNR = $snr, local")
end
labels = permutedims(labels)
plot!(
    1:20,
    p,
    xlabel="Iteration",
    label=labels,
    lw=2,
    ls=:dash,
    title=titlefer,
    ylims=(lim,0)
)

p = readdlm("./Saved Data/Clipping/FER_List-RBP 0.8 clipping.txt",'\t',Float64,'\n')
labels = Vector{String}()
for snr in SNR
    push!(labels,"SNR = $snr, clipping")
end
labels = permutedims(labels)
plot!(
    1:20,
    p,
    xlabel="Iteration",
    label=labels,
    lw=2,
    ls=:dashdot,
    title=titlefer,
    ylims=(lim,0)
)

