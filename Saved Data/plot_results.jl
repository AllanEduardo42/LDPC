p1 = readdlm("./Saved Data/2025-02-18T17:08:56.641/FER_List-RBP 0.8.txt",'\t',Float64,'\n')

p2 = readdlm("./Saved Data/2025-02-19T17:08:59.761/FER_List-RBP 0.8.txt",'\t',Float64,'\n')

p = readdlm("./Saved Data/2025-02-18T17:08:56.641/FER_RBP 0.8.txt",'\t',Float64,'\n')

plotlyjs()
lim = log10(1/maximum(TRIALS))
labels = Vector{String}()
for snr in SNR
    push!(labels,"SNR (dB) = $snr, listsize = 2")
end
labels = permutedims(labels)
titlefer = "FER (decay factor 0.8)"
plot(
    1:16,
    p1,
    xlabel="Iteration",
    label=labels,
    lw=2,
    title=titlefer,
    ylims=(lim,0)
)
labels = Vector{String}()
for snr in SNR
    push!(labels,"SNR (dB) = $snr, listsize = 1")
end
labels = permutedims(labels)
plot!(
    1:16,
    p2,
    xlabel="Iteration",
    label=labels,
    lw=2,
    ls=:dash,
    title=titlefer,
    ylims=(lim,0)
)
labels = Vector{String}()
for snr in SNR
    push!(labels,"SNR (dB) = $snr, no list")
end
labels = permutedims(labels)
plot!(
    1:16,
    p,
    xlabel="Iteration",
    label=labels,
    lw=2,
    ls=:dot,
    title=titlefer,
    ylims=(lim,0)
)