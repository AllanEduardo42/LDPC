using DelimitedFiles
using Plots

N = 576
R = [1,2]
EbNo = 2.25
protocol = "NR5G"
decays = collect(0.1:0.1:1.0)
FB = ["F","B"]
maxiter = 50

algorithms = ["SVNF",
               "D-SVNF"]

plotlyjs()
directory = "./Saved Data/D-SVNF $protocol $N $(R[1])|$(R[2]) $EbNo/"
lim = log10(1/maximum(12800000))+1

for j=1:2
    title = FB[j]*"ER $protocol (N = $N, R = $(R[1])/$(R[2]), Eb/N0 = $(EbNo)dB)"
    p = plot()
    for i in eachindex(algorithms)
        x = readdlm(directory*FB[j]*"ER_"*algorithms[i]*".txt",'\t',Float64,'\n')
        x = x[1:maxiter,:]
        labels = Vector{String}()
        push!(labels,algorithms[i])
        labels = permutedims(labels)
        p = plot!(
            1:maxiter,
            log10.(x),
            xlabel="Iteration",
            label=labels,
            lw=2,
            title=title,
            ylims=(lim,0),
            xlim=(1,maxiter),
            color=i,
            legend_title = "Algotithm",
            legend_title_font_pointsize = 8,
            # legend = :outertopright,
            size = 1.2*[600,400]
        )

    end
    display(p)
    # Plots.pdf(p,"rbp.pdf")
    global lim -= 1
end

