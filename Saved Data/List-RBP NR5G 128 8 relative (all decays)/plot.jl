using DelimitedFiles
using Plots

SNR = [1.2, 1,4, 1.6, 1.8]
decays = [0.7, 0.8, 0.9, 1.0]
# decays = collect(0.0:0.1:1.0)

lines = [:dash,:solid,:dashdot,:dot,:dash,:solid,:dashdot,:dot,:dash,:solid,:dashdot,:dot]
FB = ["F","B"]
maxiter = 10

plotlyjs()
directory = "./Saved Data/List-RBP NR5G 128 8 relative (all decays)/"
lim = log10(1/maximum(1000*2^10))+1

for j=1:2
    title = FB[j]*"ER List-RBP (128,8) Relative"
    p = plot()
    for i in eachindex(decays)
        sdecay = string(decays[i])
        x = readdlm(directory*FB[j]*"ER_List-RBP "*sdecay*".txt",'\t',Float64,'\n')
        x = x[1:maxiter,:]
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
            xlim=(1,maxiter),
            color=[1 2 3 4],
            legend_title = "(decay, SNR)",
            legend_title_font_pointsize = 8,
            legend = :outertopright,
            size = (900,600)
        )
    end
    display(p)
    # Plots.pdf(p,"rbp.pdf")
    global lim -= 1
end

