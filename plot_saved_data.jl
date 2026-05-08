using DelimitedFiles
using Plots
using LaTeXStrings

# Date = "2026-05-05"
# Hour = "08:50"
# Protocol = "5GNR"
# N = 576
# R = [1,2]
# EbN0 = [1.0, 1.5, 2.0, 2.5, 3.0]

Date = "2026-05-04"
Hour = "11:26"
Protocol = "5GNR"
N = 1248
R = [1,2]
EbN0 = [1.0, 1.25, 1.5, 1.75, 2.0]

maxiter = 20
iter = 10

FB = ["F","B"]
# markers = [:none, :none, :dtriangle, :circle, :rect, :utriangle, :diamond, :cross, :star5, :hexagon]
modes_markers_lines =  ["Flooding"             :none            :solid          
                        "LBP"                  :none            :solid
                        "RBP"                  :dtriangle       :solid
                        "RD-RBP"               :utriangle       :solid                
                        "NW-RBP"               :star5           :solid         
                        "SVNF"                 :diamond         :solid
                        "List-RBP"             :circle          :solid
                        "UBP-RBP"              :cross            :solid
                        "C-RBP"                :none            :dot
                        "C&R-RBP"              :rect            :solid
                        "C&DR-RBP"             :circle          :solid
                        # "C&DR-RBP 4"             :circle          :solid
                 ]


directory = "./Saved Data/$Date $Hour $Protocol $N $(R[1])|$(R[2])"
for ebn0 in EbN0
    global directory *= " $(ebn0)dB"
end
directory *= "/"

FER_EbN0 = zeros(length(EbN0),size(modes_markers_lines,1))
BER_EbN0 = zeros(length(EbN0),size(modes_markers_lines,1))

p1 = plot()
p2 = plot()


colors = [1, 2, 3, 4, 5, 15, :gray, :black, 11, 11, 14, 13, 14]
marker_sizes = [4,4,4,4,5,5,3,3,3,4,3,4,4]


# #
gr()

values = 10.0 .^(-[0, 1, 2, 3, 4, 5])
labels = ["0", "-1", "-2", "-3", "-4", "-5", "-6"]

min_x = zeros(size(modes_markers_lines,1))
max_x = zeros(size(modes_markers_lines,1))

for k in eachindex(EbN0)
    for j=1:2
        title = FB[j]*"ER $Protocol (N = $N, R = $(R[1])/$(R[2]), Eb/N0 = $(EbN0[k])dB)"
        p = plot()
        for i in axes(modes_markers_lines,1)
            str = modes_markers_lines[i,1]
            x = readdlm(directory*FB[j]*"ER_"*str*".txt",'\t',Float64,'\n')
            x = x[1:maxiter,k]
            min_x[i] = minimum(x)
            max_x[i] = maximum(x)
            p = plot!(
                1:maxiter,
                # log10.(x),
                x,
                ylabel=FB[j]*"ER",
                label=modes_markers_lines[i,1],
                lw=1.5,
                ls=modes_markers_lines[i,3],
                # title=title,
                # ylims=(liminf,limsup),
                # xlim=(1,maxiter),
                minorgrid=true,
                minorgridalpha=0.05,
                yscale=:log10,
                color=colors[i],
                markersize=marker_sizes[i],
                markerstrokewidth = 0.5,
                guidefontsize=10,
                # ticks=9,
                # yticks = (values,labels),
                # tickfontsize=10,
                # legend_title = "(Algorithm, best decay factor)",
                # legend_title_font_pointsize = 10,
                # legend_font_pointsize = 10,
                legend_position = :topright,
                # legend = :outertopright,
                size = 1.5 .*(300,500),
                markershape = modes_markers_lines[i,2],
                left_margin=3Plots.mm,
                # bottom_margin=-20Plots.mm,
                top_margin=1Plots.mm,
                fontfamily="Computer Modern",
                framestyle=:box
            )
        end
        liminf = 10^(floor(log10(minimum(min_x))))
        limsup = 10^(ceil(log10(maximum(max_x))))
        p = plot!(ylims=(liminf,limsup))
        if j == 1
            global p1 = p
        else
            p2 = plot!(xlabel="Iteration")
            global p2 = p
        end
    end
    
    p = plot(p1,p2,layout=(2,1),bottom_margin=-3Plots.mm)
    Plots.pdf(p,directory*"/FER_BER_$(N)_$(EbN0[k])dB_$maxiter")
end

# FER x EbN0
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
    ylims=(10^(-4),1),
    color=permutedims(colors),
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
    size = 1.5 .*(300,300),
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
Plots.pdf(PLOT,directory*"/FER_$(iter)_dB.pdf")

# #

# plotlyjs()

# for k in eachindex(EbN0)
#     for j=1:2
#         title = FB[j]*"ER $Protocol (N = $N, R = $(R[1])/$(R[2]), Eb/N0 = $(EbN0[k])dB)"
#         p = plot()
#         for i in axes(modes_markers_lines,1)
#             str = modes_markers_lines[i,1]
#             x = readdlm(directory*FB[j]*"ER_"*str*".txt",'\t',Float64,'\n')
#             x = x[:,k]
#             p = plot!(
#                 1:maxiter,
#                 log10.(x[1:maxiter]),
#                 lw=4,
#                 title=title,
#                 label=modes_markers_lines[i,1],
#                 markershape = modes_markers_lines[i,2],
#                 size = (1400,1000),
#                 color=colors[i],
#                 tickfontsize=16,
#                 legend_font_pointsize = 20,
#                 markersize=8,
#                 markerstrokewidth = 1
#             )
#         end
#         display(p)
#     end    
# end

# # #
# # FER x EbN0
# if length(EbN0) > 1
#     for j = 1:2
#         for k in eachindex(EbN0)
#             for i in axes(modes_markers_lines,1)
#                 str = modes_markers_lines[i,1]
#                 x = readdlm(directory*FB[j]*"ER_"*str*".txt",'\t',Float64,'\n')
#                 x = x[:,k]
#                 if j == 1
#                     FER_EbN0[k,i] = log10.(x[iter])
#                 else
#                     BER_EbN0[k,i] = log10.(x[iter])
#                 end
#             end
#         end
#     end

#     p = plot(
#         EbN0,FER_EbN0,
#         title = "FER x EbN0 $Protocol (N = $N, R = $(R[1])/$(R[2]), Iter = $iter)",
#         color=permutedims(colors),
#         label=permutedims(modes_markers_lines[:,1]),
#         markershape=permutedims(modes_markers_lines[:,2]),
#         line=permutedims(modes_markers_lines[:,3]),
#         lw=4,
#         markersize=8,
#         markerstrokewidth = 1,
#         guidefontsize=10,
#         legend_position = :bottomleft,
#         tickfontsize=16,
#         legend_font_pointsize = 20,
#         size = (1400,1000)
#     )
#     display(p)
# end
