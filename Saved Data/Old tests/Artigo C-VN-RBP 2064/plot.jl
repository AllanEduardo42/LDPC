N = 2064
R = [1,2]
iters = 1024000
maxiter = 50
EbNo = 1.5
protocol = "NR5G"

using DelimitedFiles
using Plots
using LaTeXStrings

function save_pdf(p, filename)
    width, height = Int.(p.attr[:size])
    Plots.prepare_output(p)
    PlotlyJS.savefig(Plots.plotlyjs_syncplot(p), filename*".pdf", width=width, height=height)
end

plotlyjs()
# gr()

FB = ["F","B"]
# markers = [:none, :none, :dtriangle, :circle, :rect, :utriangle, :diamond, :cross, :star5, :hexagon]
modes_markers = ["VN-RBP 0.85"          :circle 
                 "C-VN-RBP 0.85 3+"     :rect
                 "C-VN-RBP 0.85 3"      :dtriangle
                 ]
                 
directory = "./Saved Data/Artigo C-VN-RBP 2064/"
liminf = 10^(-3)
limsup = 1

for j=1:2
    title = FB[j]*"ER $protocol (N = $N, R = $(R[1])/$(R[2]), Eb/N0 = $(EbNo)dB)"
    p = plot()
    for k in axes(modes_markers,1)
        str = modes_markers[k,1]
        labels = modes_markers[k,1]
        x = readdlm(directory*FB[j]*"ER_"*str*".txt",'\t',Float64,'\n')
        x = x[1:maxiter,:]
        line = :solid
        # labels = permutedims(labels)
        p = plot!(
            1:maxiter,
            x,
            xlabel="Iteration",
            ylabel=FB[j]*"ER",
            label=labels,
            lw=3,
            ls=line,
            title=title,
            ylims=(liminf,limsup),
            # xlim=(1,maxiter),
            minorgrid=true,
            yscale=:log10,
            color=k,
            markersize=5,
            guidefontsize=20,
            tickfontsize=15,
            # legend_title = "(Algorithm, best decay factor)",
            # legend_title_font_pointsize = 10,
            legend_font_pointsize = 15,
            # legend = :outertopright,
            size = 1.5 .*(600,400),
            markershape = modes_markers[k,2],
            left_margin=3Plots.mm,
            bottom_margin=3Plots.mm,
            top_margin=3Plots.mm
        )
    end
    display(p)
    # save_pdf(p,directory*"/$(FB[j])ER")
    # Plots.pdf(p,directory*"/$(FB[j])ER")
    global liminf = 10^(-4)
    global limsup = 10^(-2)
end







