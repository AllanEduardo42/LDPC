N = 576
R = [1,2]
iters = 1_280_000
maxiter = 20
iter = 5
EbN0 = [2.5]
active = [1]
protocol = "NR5G"
decays = [0.7, 0.8, 0.9, 1.0, 0.85]

plot_mode = "gr"
# plot_mode = "plotlyjs"

using DelimitedFiles
using Plots
using LaTeXStrings

if plot_mode == "gr"
    gr()
else
    plotlyjs()
end

FB = ["F","B"]
# markers = [:none, :none, :dtriangle, :circle, :rect, :utriangle, :diamond, :cross, :star5, :hexagon]
modes =  ["RD-RBP", "List-RBP (16,2)", "C&R-RBP", "C&DR-RBP", "C-RBP"]
                 
directory = "./Saved Data/Decays response letter (5G)/"

markers = [:diamond, :circle, :rect, :cross, :utriangle ]

marker_sizes = [4,4,4,5,4]


for i in axes(modes,1)
    p1 = plot()
    p2 = plot()
    if plot_mode == "gr"
        liminf = 1*10^(-4)
        limsup = 1*1*10^(-0)
    else
        liminf = -4
        limsup = -0
    end

    for k = 1:2        
        p = plot()
        for j in eachindex(decays)
            str = modes[i,1]*" $(decays[j])"
            x = readdlm(directory*FB[k]*"ER_"*str*".txt",'\t',Float64,'\n')
            if j == 5
                x = x[:,4] # column that refers to 2.5 dB in the file
            end
            if plot_mode != "gr"
                x = log10.(x)
            end
            p = plot!(
                1:maxiter,
                # log10.(x),
                x[1:maxiter],
                ylabel=FB[k]*"ER",
                label=str,
                lw=2,
                # ls=modes_markers_lines[i,3],
                # title=title,
                ylims=(liminf,limsup),
                # xlim=(1,maxiter),
                minorgrid=true,
                minorgridalpha=0.05,
                yscale=:log10,
                markersize=marker_sizes[j],
                markerstrokewidth = 0.5,
                guidefontsize=10,
                # tickfontsize=10,
                # legend_title = "(Algorithm, best decay factor)",
                # legend_title_font_pointsize = 10,
                # legend_font_pointsize = 10,
                legend_position = :topright,
                # legend = :outertopright,
                size = 1.5 .*(300,400),
                markershape = markers[j],
                left_margin=1Plots.mm,
                # bottom_margin=-20Plots.mm,
                top_margin=1Plots.mm,
                fontfamily="Computer Modern",
                framestyle=:box
            )
        end
        if k == 1
            p1 = p
        else
            p2 = plot!(xlabel="Iteration")
            p2 = p
        end
        if plot_mode == "gr"
            liminf = 1*10^(-5)
            limsup = 1*10^(-1)
        else
            liminf = -5
            limsup = -1
        end
    end
    p = plot(p1,p2,layout=(2,1),bottom_margin=-3Plots.mm,)
    # display(p)
    Plots.pdf(p,directory*"/FER_BER_576_$(modes[i])")
end
