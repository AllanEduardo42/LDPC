N = 2064
R = [1,2]
iters = 512000
maxiter = 50
EbNo = 1.5
protocol = "NR5G"

using DelimitedFiles
using Plots

function save_pdf(p, filename)
    width, height = Int.(p.attr[:size])
    Plots.prepare_output(p)
    PlotlyJS.savefig(Plots.plotlyjs_syncplot(p), filename*".pdf", width=width, height=height)
end

plotlyjs()

FB = ["F","B"]
markers = [:none, :none, :dtriangle, :circle, :rect, :utriangle, :diamond, :cross,:star5]
modes = ["RBP","E-RBP","FB-RBP","VN-RBP"]
directory = "./Saved Data/Testes EFB $protocol $N $(R[1])|$(R[2]) $EbNo/"
liminf = -3.5
limsup = 0.1

for j=1:2
    title = FB[j]*"ER $protocol (N = $N, R = $(R[1])/$(R[2]), Eb/N0 = $(EbNo)dB)"
    p = plot()
    for k in eachindex(modes)
        str = modes[k]*" 1.0"
        labels = modes[k]
        x = readdlm(directory*FB[j]*"ER_"*str*".txt",'\t',Float64,'\n')
        x = x[1:maxiter,:]
        line = :solid
        # labels = permutedims(labels)
        p = plot!(
            1:maxiter,
            log10.(x),
            xlabel="Iteration",
            ylabel="log (FER)",
            label=labels,
            lw=2,
            ls=line,
            title=title,
            ylims=(liminf,limsup),
            # xlim=(1,maxiter),
            color=k,
            legend_title = "(algorithm, best decay factor)",
            legend_title_font_pointsize = 10,
            legend_font_pointsize = 10,
            # legend = :outertopright,
            size = 1.6 .*(600,400),
            markershape=markers[k]
        )
    end
    display(p)
    save_pdf(p,directory*"/$(FB[j])ER")
    global liminf -= 1
    global limsup -= 1
end







