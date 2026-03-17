N = 576
R = [1,2]
iters = 12800000
maxiter = 20
iter = 5
EbN0 = [1.0, 1.5, 2.0, 2.5, 3.0]
active_EbN0 = [0 0 0 1 0]
protocol = "NR5G"

using DelimitedFiles
using Plots
using LaTeXStrings

function reorganize_plot_data(plot_data)

    new_plot_data = []
    for i in axes(plot_data,1)
        if plot_data[i,4]
            push!(new_plot_data,plot_data[i,:])
        end
    end

    new_plot_data = reduce(hcat,new_plot_data)

    return permutedims(new_plot_data)

end

# plotlyjs()
gr()
# pythonplot()

FB = ["F","B"]
# markers = [:none, :none, :dtriangle, :circle, :rect, :utriangle, :diamond, :cross, :star5, :hexagon]
plot_data =  
[   #algorithm             #marker      #line   #plot?  #color  #marker size
    "Flooding"             :none        :solid  false   1       4      
    "LBP"                  :none        :solid  false   2       4
    "RBP"                  :dtriangle   :solid  false   3       4
    "RD-RBP 0.85"          :utriangle   :solid  true    4       4            
    "NW-RBP"               :star5       :solid  false   5       5  
    "SVNF"                 :diamond     :solid  true    6       5
    "List-RBP (16,2) 0.85" :circle      :dot    false   7       3
    "C&R-RBP 0.85"         :rect        :solid  true    8       4
    "C&DR-RBP 0.85"        :circle      :solid  true    9       3
#    "C-RBP 0.85"           :none        :dot    false   9       3
#    "R-RBP 0.85"           :xcross      :dot    false   9       3
    "CI-RBP"                :xcross     :solid  true    2      4
]

plot_data = reorganize_plot_data(plot_data)


directory = "./Saved Data/CI-RBP 576 1|2 1:3dB/"
# liminf = -5
# limsup = 0

FER_EbN0 = zeros(length(EbN0),size(plot_data,1))
BER_EbN0 = zeros(length(EbN0),size(plot_data,1))

p1 = plot()
p2 = plot()

liminf = 1*10^(-4)
limsup = 1

for j=1:2
    for k in eachindex(EbN0)
        if active_EbN0[k] == 1
            title = FB[j]*"ER $protocol (N = $N, R = $(R[1])/$(R[2]), Eb/N0 = $(EbN0[k])dB)"
            p = plot()
            for i in axes(plot_data,1)
                str = plot_data[i,1]
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
                    label=plot_data[i,1],
                    lw=2,
                    ls=plot_data[i,3],
                    # title=title,
                    ylims=(liminf,limsup),
                    # xlim=(1,maxiter),
                    minorgrid=true,
                    minorgridalpha=0.05,
                    yscale=:log10,
                    color=plot_data[i,5],
                    markersize=plot_data[i,6],
                    markerstrokewidth = 0.5,
                    guidefontsize=10,
                    # tickfontsize=10,
                    # legend_title = "(Algorithm, best decay factor)",
                    # legend_title_font_pointsize = 10,
                    # legend_font_pointsize = 10,
                    legend_position = :topright,
                    # legend = :outertopright,
                    size = 1.5 .*(300,400),
                    markershape = plot_data[i,2],
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
        for i in axes(plot_data,1)
            str = plot_data[i,1]
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
    color=plot_data[:,5]',
    label=permutedims(plot_data[:,1]),
    markershape=permutedims(plot_data[:,2]),
    line=permutedims(plot_data[:,3]),
    lw=1.5,
    markersize=plot_data[:,6]',
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
# BER x EbN0
# PLOT = plot(
#     EbN0,BER_EbN0,
#     minorgrid=true,
#     yscale=:log10,
#     xlabel="EbN0 (dB)",
#     ylims=(10^(-8),0.1),
#     label=permutedims(plot_data[:,1]),
#     markershape=permutedims(plot_data[:,2]),
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

for k in axes(plot_data,1)
    if plot_data[k,4]
        str = plot_data[k,1]
        labels = plot_data[k,1]
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
            ls=plot_data[k,3],
            # title=title,
            ylims=(liminf,limsup),
            # xlim=(1,maxiter),
            minorgrid=true,
            minorgridalpha=0.05,
            yscale=:log10,
            color=plot_data[k,5],
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
            markershape = plot_data[k,2],
            left_margin=5Plots.mm,
            bottom_margin=5Plots.mm,
            top_margin=3Plots.mm,
            framestyle=:box
        )
    end
end

# p = plot(p1,p2,layout=(2,1),fontfamily="Computer Modern",xlabel="Iteration",bottom_margin=-1Plots.mm,)
# display(p)
# Plots.pdf(p,directory*"/FER_BER_C&DR_576")