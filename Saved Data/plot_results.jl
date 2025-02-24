using DelimitedFiles
using Plots

SNR = collect(1.2:0.4:2.0)

p1 = readdlm("./Saved Data/2025-02-20T12:25:21.943/FER_List-RBP 0.8.txt",'\t',Float64,'\n')

p2 = readdlm("./Saved Data/2025-02-20T12:25:21.943/FER_RBP 0.8.txt",'\t',Float64,'\n')

p3 = readdlm("./Saved Data/2025-02-20T12:25:21.943/FER_Local-RBP 0.8.txt",'\t',Float64,'\n')

p4 = readdlm("./Saved Data/2025-02-21T11:01:04.780/FER_List-RBP 0.8.txt",'\t',Float64,'\n')



plotlyjs()
lim = log10(1/maximum(1000*2^10))
labels = Vector{String}()
for snr in SNR
    push!(labels,"SNR (dB) = $snr, listsize = 1")
end
labels = permutedims(labels)
titlefer = "FER (decay factor 0.8)"
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
    push!(labels,"SNR (dB) = $snr, no list")
end
labels = permutedims(labels)
plot!(
    1:20,
    p2,
    xlabel="Iteration",
    label=labels,
    lw=2,
    ls=:dot,
    title=titlefer,
    ylims=(lim,0)
)

labels = Vector{String}()
for snr in SNR
    push!(labels,"SNR (dB) = $snr, local")
end
labels = permutedims(labels)
plot!(
    1:20,
    p3,
    xlabel="Iteration",
    label=labels,
    lw=2,
    ls=:dash,
    title=titlefer,
    ylims=(lim,0)
)

labels = Vector{String}()
for snr in SNR
    push!(labels,"SNR (dB) = $snr, new")
end
labels = permutedims(labels)
display(plot!(
    1:20,
    p4,
    xlabel="Iteration",
    label=labels,
    lw=2,
    ls=:dashdot,
    title=titlefer,
    ylims=(lim,0)
))

labels = Vector{String}()
for snr in SNR
    push!(labels,"SNR (dB) = $snr, old")
end
labels = permutedims(labels)
plot(
    1:20,
    p2,
    xlabel="Iteration",
    label=labels,
    lw=2,
    ls=:dash,
    title=titlefer,
    ylims=(lim,0)
)
p5 = readdlm("./Saved Data/2025-02-23T12:03:26.851/FER_RBP 0.8.txt",'\t',Float64,'\n')
labels = Vector{String}()
for snr in SNR
    push!(labels,"SNR (dB) = $snr, new")
end
labels = permutedims(labels)
plot!(
    1:20,
    p5,
    xlabel="Iteration",
    label=labels,
    lw=2,
    ls=:dot,
    title=titlefer,
    ylims=(lim,0)
)
p6 = readdlm("./Saved Data/Testing update of Lq[vnmax,cnmax] in RBP/FER_RBP 0.8.txt",'\t',Float64,'\n')
labels = Vector{String}()
for snr in SNR
    push!(labels,"SNR (dB) = $snr, Lq")
end
labels = permutedims(labels)
plot!(
    1:20,
    p6,
    xlabel="Iteration",
    label=labels,
    lw=2,
    title=titlefer,
    ylims=(lim,0)
)

