N = 576
R = [1,2]
iters = 1024000
maxiter = 10
EbN0 = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5]
protocol = "NR5G"

using DelimitedFiles
using Plots
using LaTeXStrings

function save_pdf(p, filename)
    width, height = Int.(p.attr[:size])
    Plots.prepare_output(p)
    PlotlyJS.savefig(Plots.plotlyjs_syncplot(p), filename*".pdf", width=width, height=height)
end

# plotlyjs()
gr()

FB = ["F","B"]
# markers = [:none, :none, :dtriangle, :circle, :rect, :utriangle, :diamond, :cross, :star5, :hexagon]
# modes_markers = [
#                 "Flooding"             :none            
#                  "LBP"                  :none
#                  "RD-RBP 0.85"             :dtriangle
#                  "NW-RBP"               :circle
#                  "SVNF"                 :rect 
#                  "List-RBP 0.85"        :utriangle
#                  "C&R-RBP 0.85"          :diamond 
#                  "C-RBP 0.85"           :star5
#                 #  "C&DR-RBP 0.85"        :cross 
#                  ]
modes_markers = ["Flooding"             :none            
                 "LBP"                  :none
                #  "RBP"                  :dtriangle
                 "RD-RBP 0.85"          :utriangle                 
                 "NW-RBP"               :diamond
                 "SVNF"                 :circle
                 "List-RBP (16,2) 0.85" :cross
                #  "List-RBP (16,2) 1.0"  :utriangle
                 "C&R-RBP 0.85"          :rect
                #  "C-VN-RBP 0.85 no-opt"        :star5
                 ]
directory = "./Saved Data/Artigo EbN0 576/"
# liminf = -5
# limsup = 0

FERMAX = zeros(length(EbN0),size(modes_markers,1))
BERMAX = zeros(length(EbN0),size(modes_markers,1))

for j=1:2
    for ebn0 in eachindex(EbN0)
        title = FB[j]*"ER $protocol (N = $N, R = $(R[1])/$(R[2]), Eb/N0 = $(EbN0[ebn0])dB)"
        p = plot()
        for k in axes(modes_markers,1)
            str = modes_markers[k,1]
            labels = modes_markers[k,1]
            x = readdlm(directory*FB[j]*"ER_"*str*".txt",'\t',Float64,'\n')
            x = x[1:maxiter,ebn0]
            if j == 1
                FERMAX[ebn0,k] = x[maxiter]
            else
                BERMAX[ebn0,k] = x[maxiter]
            end
            line = :solid
            # labels = permutedims(labels)
            # p = plot!(
            #     1:maxiter,
            #     # log10.(x),
            #     x,
            #     xlabel="Iteration",
            #     ylabel=FB[j]*"ER",
            #     label=labels,
            #     lw=3,
            #     ls=line,
            #     # title=title,
            #     # ylims=(liminf,limsup),
            #     # xlim=(1,maxiter),
            #     minorgrid=true,
            #     yscale=:log10,
            #     color=k,
            #     markersize=5,
            #     guidefontsize=20,
            #     tickfontsize=15,
            #     # legend_title = "(Algorithm, best decay factor)",
            #     # legend_title_font_pointsize = 10,
            #     legend_font_pointsize = 15,
            #     # legend = :outertopright,
            #     size = 1.5 .*(600,400),
            #     markershape = modes_markers[k,2],
            #     left_margin=3Plots.mm,
            #     bottom_margin=3Plots.mm,
            #     top_margin=3Plots.mm
            # )
        end
        # display(p)
        # save_pdf(p,directory*"/$(FB[j])ER")
        # Plots.pdf(p,directory*"/$(FB[j])ER_"*"$(EbN0[ebn0])dB.pdf")
    end
end

colors = [1, 2, 4, 5, 6, 18, 7]

# FER x EbN0
PLOT = plot(
    EbN0,FERMAX,
    minorgrid=true,
    minorgridalpha=0.05,
    yscale=:log10,
    ylabel="FER",
    xlabel="EbN0 (dB)",
    ylims=(10^(-6),1),
    color=colors',
    label=permutedims(modes_markers[:,1]),
    markershape=permutedims(modes_markers[:,2]),
    lw=2.5,
    markersize=5,
    markerstrokewidth = 0.25,
    guidefontsize=15,
    tickfontsize=12,
    legend_font_pointsize = 10,
    fontfamily="Computer Modern",
    size = 1.5 .*(600,400),
    left_margin=3Plots.mm,
    framestyle=:box
)
display(PLOT)
# BER x EbN0
# PLOT = plot(
#     EbN0,BERMAX,
#     minorgrid=true,
#     yscale=:log10,
#     xlabel="EbN0 (dB)",
#     ylims=(10^(-8),0.1),
#     label=permutedims(modes_markers[:,1]),
#     markershape=permutedims(modes_markers[:,2]),
#     lw=2,
#     # title="BER x EbN0 (Iter = 10, N = $N, R = $(R[1])/$(R[2]))",
#     size = 1.5 .*(600,400),
#     # ylims=(-6,-1)
# )
# display(PLOT)
Plots.pdf(PLOT,directory*"/FER_dB.pdf")









