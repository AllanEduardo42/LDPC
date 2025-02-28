using DelimitedFiles
using Plots

SNR = [1.2, 1.6, 1.8, 2.0]

plotlyjs()
lim = log10(1/maximum(1000*2^10))

p = readdlm("./Saved Data/RBP decay testing/FER_RBP 1.0.txt",'\t',Float64,'\n')
labels = Vector{String}()
for snr in SNR
    push!(labels,"SNR (dB) = $snr, 1.0")
end
labels = permutedims(labels)
titlefer = "FER RBP"
plot(
    1:20,
    p,
    xlabel="Iteration",
    label=labels,
    lw=2,
    title=titlefer,
    ylims=(lim,0)
)

p = readdlm("./Saved Data/RBP decay testing/FER_RBP 0.9.txt",'\t',Float64,'\n')
labels = Vector{String}()
for snr in SNR
    push!(labels,"SNR (dB) = $snr, 0.9")
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
    ylims=(lim,0),
    color=[1 2 3 4]
)
p = readdlm("./Saved Data/RBP decay testing/FER_RBP 0.8.txt",'\t',Float64,'\n')
labels = Vector{String}()
for snr in SNR
    push!(labels,"SNR (dB) = $snr, 0.8")
end
labels = permutedims(labels)
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

# BER
lim -= 1
p = readdlm("./Saved Data/RBP decay testing/BER_RBP 1.0.txt",'\t',Float64,'\n')
labels = Vector{String}()
for snr in SNR
    push!(labels,"SNR (dB) = $snr, 1.0")
end
labels = permutedims(labels)
titlefer = "BER RBP"
plot(
    1:20,
    p,
    xlabel="Iteration",
    label=labels,
    lw=2,
    title=titlefer,
    ylims=(lim,0)
)

p = readdlm("./Saved Data/RBP decay testing/BER_RBP 0.9.txt",'\t',Float64,'\n')
labels = Vector{String}()
for snr in SNR
    push!(labels,"SNR (dB) = $snr, 0.9")
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
    ylims=(lim,0),
    color=[1 2 3 4]
)

p = readdlm("./Saved Data/RBP decay testing/BER_RBP 0.8.txt",'\t',Float64,'\n')
labels = Vector{String}()
for snr in SNR
    push!(labels,"SNR (dB) = $snr, 0.9")
end
labels = permutedims(labels)
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

