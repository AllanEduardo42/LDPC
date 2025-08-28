N = 576
R = [1,2]
iters = 1024000
maxiter = 30
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
modes_markers = ["Flooding"             :none            
                 "LBP"                  :none
                 "RBP"                  :dtriangle
                 "RD-RBP 0.85"          :utriangle                 
                 "NW-RBP"               :diamond
                 "SVNF"                 :circle
                 "List-RBP (16,2) 0.85" :cross
                 "C&R-RBP 0.85"         :rect
                 "C-RBP 0.85"           :star5
                 "C&DR-RBP 0.85"        :star5
                 ]
directory = "./Saved Data/Artigo 576/"
liminf = 10^(-4)
limsup = 1.2

colors = [1, 2, 3, 4, 5, 6, 18, 7, 8, 9]

p1 = plot()
p2 = plot()

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
            ylabel=FB[j]*"ER",
            label=labels,
            lw=2,
            ls=line,
            # title=title,
            ylims=(liminf,limsup),
            # xlim=(1,maxiter),
            minorgrid=true,
            minorgridalpha=0.05,
            yscale=:log10,
            color=colors[k],
            markersize=5,
            markerstrokewidth = 0.25,
            guidefontsize=15,
            tickfontsize=12,
            # legend_title = "(Algorithm, best decay factor)",
            # legend_title_font_pointsize = 10,
            legend_font_pointsize = 10,
            # legend = :outertopright,
            size = 1.5 .*(600,500),
            markershape = modes_markers[k,2],
            left_margin=3Plots.mm,
            bottom_margin=-3Plots.mm,
            top_margin=3Plots.mm,
            framestyle=:box
        )
    end
    if j == 1
        global p1 = p
    else
        global p2 = p
    end
    # display(p)
    # # save_pdf(p,directory*"/$(FB[j])ER")
    global liminf = 10^(-5)
    global limsup = 0.1
end

# display(p1)
# display(p2)
p = plot(p1,p2,layout=(2,1),fontfamily="Computer Modern",xlabel="Iteration",bottom_margin=-1Plots.mm,)
display(p)
Plots.pdf(p,directory*"/FER_BER_576")

### C&DR-RBP

modes_markers = [           
                 "C&R-RBP 0.85"         :rect
                #  "C-RBP 0.85"           :circle
                 "C&DR-RBP 0.85"        :circle
                 ]

p = plot()
liminf = 10^(-4)
limsup = 1.2

for k in axes(modes_markers,1)
    str = modes_markers[k,1]
    labels = modes_markers[k,1]
    x = readdlm(directory*"FER_"*str*".txt",'\t',Float64,'\n')
    x = x[1:maxiter,:]
    line = :solid
    # labels = permutedims(labels)
    p = plot!(
        1:maxiter,
        x,
        ylabel="FER",
        label=labels,
        lw=2,
        ls=line,
        # title=title,
        ylims=(liminf,limsup),
        # xlim=(1,maxiter),
        minorgrid=true,
        minorgridalpha=0.05,
        yscale=:log10,
        color=colors[k],
        markersize=5,
        markerstrokewidth = 0.25,
        guidefontsize=15,
        tickfontsize=12,
        fontfamily="Computer Modern",
        xlabel="Iteration",
        # legend_title = "(Algorithm, best decay factor)",
        # legend_title_font_pointsize = 10,
        legend_font_pointsize = 10,
        # legend = :outertopright,
        size = 1.5 .*(600,300),
        markershape = modes_markers[k,2],
        left_margin=5Plots.mm,
        bottom_margin=5Plots.mm,
        top_margin=3Plots.mm,
        framestyle=:box
    )
end

# p = plot(p1,p2,layout=(2,1),fontfamily="Computer Modern",xlabel="Iteration",bottom_margin=-1Plots.mm,)
display(p)
Plots.pdf(p,directory*"/FER_BER_C&DR_576")





