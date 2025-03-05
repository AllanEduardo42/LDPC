using DelimitedFiles
using Plots

SNR = [1.2, 1.6, 1.8, 2.0]
decays = [0.7, 0.8, 0.9, 1.0]
lines = [:dashdot, :dot, :dash, :solid]
FB = ["F","B"]
maxiter = 30

plotlyjs()
directory = "./Saved Data/Local-RBP decay test/"
lim = log10(1/maximum(1000*2^10))

for j=1:2
    title = FB[j]*"ER Local-RBP"
    p = plot()
    for i in eachindex(decays)
        sdecay = string(decays[i])
        x = readdlm(directory*FB[j]*"ER_Local-RBP "*sdecay*".txt",'\t',Float64,'\n')
        labels = Vector{String}()
        for snr in SNR
            push!(labels,sdecay*", $(snr)dB")
        end
        labels = permutedims(labels)
        p = plot!(
            1:maxiter,
            x,
            xlabel="Iteration",
            label=labels,
            lw=2,
            ls=lines[i],
            title=title,
            ylims=(lim,0),
            color=[1 2 3 4],
            legend_title = "(decay, SNR)",
            legend_title_font_pointsize = 8,
            legend = :outertopright,
            size = (900,600)
        )
    end
    display(p)
    global lim -= 1
end

