N = 576
R = [1,2]
iters = 12800000
maxiter = 20
iter = 10
EbN0 = [1.0, 1.5, 2.0, 2.5, 3.0]
active = [0 0 0 1 0]
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
# pythonplot()

FB = ["F","B"]
# markers = [:none, :none, :dtriangle, :circle, :rect, :utriangle, :diamond, :cross, :star5, :hexagon]
modes_markers_lines =  ["Flooding"             :none            :solid          
                        "LBP"                  :none            :solid
                        "RBP"                  :dtriangle       :solid
                        "RD-RBP 0.85"          :utriangle       :solid                
                        "NW-RBP"               :star5         :solid         
                        "SVNF"                 :diamond          :solid
                        "List-RBP (16,2) 0.85" :circle        :dot
                        "C&R-RBP 0.85"         :rect            :solid
                        "C&DR-RBP 0.85"        :circle        :solid
                        #  "C-RBP 0.85"          :none           :dot
                        #  "R-RBP 0.85"           :xcross
                 ]
directory = "./Saved Data/Testes NR5G 576 1|2 1.0:3.0dB/"
# liminf = -5
# limsup = 0

FER_EbN0 = zeros(length(EbN0),size(modes_markers_lines,1))
BER_EbN0 = zeros(length(EbN0),size(modes_markers_lines,1))

p1 = plot()
p2 = plot()

liminf = 1*10^(-4)
limsup = 1
colors = [1, 2, 3, 4, 5, 6, 7, 8, 9, 9]
marker_sizes = [4,4,4,4,5,5,3,4,3,3]

for j=1:2
    for k in eachindex(EbN0)
        if active[k] == 1
            title = FB[j]*"ER $protocol (N = $N, R = $(R[1])/$(R[2]), Eb/N0 = $(EbN0[k])dB)"
            p = plot()
            for i in axes(modes_markers_lines,1)
                str = modes_markers_lines[i,1]
                x = readdlm(directory*FB[j]*"ER_"*str*".txt",'\t',Float64,'\n')
                x = x[:,k]
                if j == 1
                    FER_EbN0[k,i] = x[iter]
                else
                    BER_EbN0[k,i] = x[iter]
                end
                p = plot!(
                    1:maxiter,
                    # log10.(x),
                    x[1:maxiter],
                    ylabel=FB[j]*"ER",
                    label=modes_markers_lines[i,1],
                    lw=2,
                    ls=modes_markers_lines[i,3],
                    # title=title,
                    ylims=(liminf,limsup),
                    # xlim=(1,maxiter),
                    minorgrid=true,
                    minorgridalpha=0.05,
                    yscale=:log10,
                    color=colors[i],
                    markersize=marker_sizes[i],
                    markerstrokewidth = 0.5,
                    guidefontsize=10,
                    # tickfontsize=10,
                    # legend_title = "(Algorithm, best decay factor)",
                    # legend_title_font_pointsize = 10,
                    # legend_font_pointsize = 10,
                    legend_position = :topright,
                    # legend = :outertopright,
                    size = 1.5 .*(300,400),
                    markershape = modes_markers_lines[i,2],
                    left_margin=1Plots.mm,
                    # bottom_margin=-20Plots.mm,
                    top_margin=1Plots.mm,
                    fontfamily="Computer Modern",
                    framestyle=:box
                )
            end
            if j == 1
                global p1 = p
            else
                p2 = plot!(xlabel="Iteration")
                global p2 = p
            end
            global liminf = 1*10^(-5)
            global limsup = 0.1
            # display(p)
        end
        # save_pdf(p,directory*"/$(FB[j])ER")
        # Plots.pdf(p,directory*"/$(FB[j])ER_"*"$(EbN0[k])dB.pdf")
    end
end
p = plot(p1,p2,layout=(2,1),bottom_margin=-3Plots.mm,)
# display(p)
Plots.pdf(p,directory*"/FER_BER_576_new")



# # FER x EbN0
for j = 1:2
    for k in eachindex(EbN0)
        for i in axes(modes_markers_lines,1)
            str = modes_markers_lines[i,1]
            x = readdlm(directory*FB[j]*"ER_"*str*".txt",'\t',Float64,'\n')
            x = x[:,k]
            if j == 1
                FER_EbN0[k,i] = x[iter]
            else
                BER_EbN0[k,i] = x[iter]
            end
        end
    end
end

PLOT = plot(
    EbN0,FER_EbN0,
    minorgrid=true,
    minorgridalpha=0.05,
    yscale=:log10,
    ylabel="FER",
    xlabel=L"E_b/N_0"*" (dB)",
    ylims=(10^(-5),1),
    color=colors',
    label=permutedims(modes_markers_lines[:,1]),
    markershape=permutedims(modes_markers_lines[:,2]),
    line=permutedims(modes_markers_lines[:,3]),
    lw=1.5,
    markersize=marker_sizes',
    markerstrokewidth = 0.5,
    guidefontsize=10,
    legend_position = :bottomleft,
    # tickfontsize=12,
    # legend_font_pointsize = 10,
    fontfamily="Computer Modern",
    size = 1.5 .*(300,200),
    # left_margin=1Plots.mm,
    framestyle=:box
)
# display(PLOT)
# BER x EbN0
# PLOT = plot(
#     EbN0,BER_EbN0,
#     minorgrid=true,
#     yscale=:log10,
#     xlabel="EbN0 (dB)",
#     ylims=(10^(-8),0.1),
#     label=permutedims(modes_markers_lines[:,1]),
#     markershape=permutedims(modes_markers_lines[:,2]),
#     lw=2,
#     # title="BER x EbN0 (Iter = 10, N = $N, R = $(R[1])/$(R[2]))",
#     size = 1.5 .*(600,400),
#     # ylims=(-6,-1)
# )
# display(PLOT)
Plots.pdf(PLOT,directory*"/FER_dB.pdf")

### C&DR-RBP

p = plot()
liminf = 10^(-3)
limsup = 1

for k in axes(modes_markers_lines,1)
    str = modes_markers_lines[k,1]
    labels = modes_markers_lines[k,1]
    x = readdlm(directory*"FER_"*str*".txt",'\t',Float64,'\n')
    x = x[1:8,4]
    line = :solid
    # labels = permutedims(labels)
    p = plot!(
        1:8,
        x,
        ylabel="FER",
        label=labels,
        lw=3,
        ls=modes_markers_lines[k,3],
        # title=title,
        ylims=(liminf,limsup),
        # xlim=(1,maxiter),
        minorgrid=true,
        minorgridalpha=0.05,
        yscale=:log10,
        color=colors[k],
        markersize=4,
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
        markershape = modes_markers_lines[k,2],
        left_margin=5Plots.mm,
        bottom_margin=5Plots.mm,
        top_margin=3Plots.mm,
        framestyle=:box
    )
end

# p = plot(p1,p2,layout=(2,1),fontfamily="Computer Modern",xlabel="Iteration",bottom_margin=-1Plots.mm,)
# display(p)
Plots.pdf(p,directory*"/FER_BER_C&DR_576")