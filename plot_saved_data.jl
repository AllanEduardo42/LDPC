using DelimitedFiles
using Plots
using LaTeXStrings

# Date = "2026-05-05"
# Hour = "08:50"
# Protocol = "5GNR"
# N = 576
# R = [1,2]
# EbN0 = [1.0, 1.5, 2.0, 2.5, 3.0]

# Date = "2026-05-04"
# Hour = "11:26"
# Protocol = "5GNR"
# N = 1248
# R = [1,2]
# EbN0 = [1.0, 1.25, 1.5, 1.75, 2.0]

# Date = "2026-05-08"
# Hour = "16:37"
# Protocol = "5GNR"
# N = 576
# R = [2,3]
# EbN0 = [2.0, 2.5, 3.0, 3.5, 4.0]

Date = "2026-05-08"
Hour = "22:11"
Protocol = "5GNR"
N = 576
R = [1,2]
EbN0 = [3.5, 4.0, 4.5, 5.0, 5.5, 6.0]

# Date = "2026-05-10"
# Hour = "18:23"
# Protocol = "5GNR"
# N = 1248
# R = [2,3]
# EbN0 = [2.0, 2.2, 2.4, 2.6, 2.8]

maxiter = 20
iter = 10

FB = ["F","B"]
# markers = [:none, :none, :dtriangle, :circle, :rect, :utriangle, :diamond, :cross, :star5, :hexagon]

 
curves =  
# algorithm             marker           line       color       marker size                        
[
 "Flooding"             :none            :solid     1           4  
 "LBP"                  :none            :solid     2           4
 "RBP"                  :dtriangle       :solid     3           4      
 "RD-RBP"               :utriangle       :solid     4           4         
 "NW-RBP"               :star5           :solid     5           5
 "SVNF"                 :diamond         :solid     15          5
 "List-RBP"             :circle          :solid     :gray       3
 "UBP-RBP"              :cross           :solid     :black      3
 "C&R-RBP"              :rect            :solid     11          4
 "C&DR-RBP"             :circle          :solid     14          4       
 "C-RBP"                :none            :dot       11          4
]


directory = "./Saved Data/$Date $Hour $Protocol $N $(R[1])|$(R[2])"
for ebn0 in EbN0
    global directory *= " $(ebn0)dB"
end
directory *= "/"

L = size(curves,1)-1
FER_EbN0 = zeros(length(EbN0),L) # C-RBP not included
BER_EbN0 = zeros(length(EbN0),L)

p1 = plot()
p2 = plot()

##
gr()

min_x = zeros(size(curves,1))
max_x = zeros(size(curves,1))

for k in eachindex(EbN0)
    for j=1:2
        title = FB[j]*"ER $Protocol (N = $N, R = $(R[1])/$(R[2]), Eb/N0 = $(EbN0[k])dB)"
        p = plot()
        for i in axes(curves,1)
            str = curves[i,1]
            x = readdlm(directory*FB[j]*"ER_"*str*".txt",'\t',Float64,'\n')
            x = x[1:maxiter,k]
            min_x[i] = minimum(x)
            max_x[i] = maximum(x)
            p = plot!(
                1:maxiter,
                x,
                ylabel=FB[j]*"ER",
                label=curves[i,1],
                lw=1.5,
                ls=curves[i,3],
                minorgrid=true,
                minorgridalpha=0.05,
                yscale=:log10,
                color=curves[i,4],
                markersize=curves[i,5],
                markerstrokewidth = 0.5,
                guidefontsize=10,
                legend_position = :topright,
                size = 1.5 .*(300,500),
                markershape = curves[i,2],
                left_margin=3Plots.mm,
                top_margin=1Plots.mm,
                fontfamily="Computer Modern",
                framestyle=:box
            )
        end
        aux = log10(minimum(min_x))
        diff = aux - floor(aux)
        if diff > 1 - 0.01
            liminf = 10^(floor(aux)+1)
        else
            liminf = 10^(floor(aux))
        end
        if j == 1
            limsup = 10^(ceil(log10(maximum(max_x))))
        else
            limsup = 0.1
        end
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
        for i in 1:L
            str = curves[i,1]
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
    lw=1.5,
    label=permutedims(curves[:,1]),
    markershape=permutedims(curves[:,2]),
    line=permutedims(curves[:,3]),
    color=permutedims(curves[:,4]),
    markersize=permutedims(curves[:,5]),
    markerstrokewidth = 0.5,
    guidefontsize=10,
    legend_position = :bottomleft,
    fontfamily="Computer Modern",
    size = 1.5 .*(300,300),
    framestyle=:box
)

Plots.pdf(PLOT,directory*"/FER_$(iter)_dB.pdf")

# #

plotlyjs()

for k in eachindex(EbN0)
    for j=1:2
        title = FB[j]*"ER $Protocol (N = $N, R = $(R[1])/$(R[2]), Eb/N0 = $(EbN0[k])dB)"
        p = plot()
        for i in axes(curves,1)
            str = curves[i,1]
            x = readdlm(directory*FB[j]*"ER_"*str*".txt",'\t',Float64,'\n')
            x = x[:,k]
            p = plot!(
                1:maxiter,
                log10.(x[1:maxiter]),
                lw=4,
                title=title,
                label=curves[i,1],
                markershape = curves[i,2],
                size = (1400,1000),
                color=curves[i,4],
                tickfontsize=16,
                legend_font_pointsize = 20,
                markersize=8,
                markerstrokewidth = 1
            )
        end
        display(p)
    end    
end

# #
# FER x EbN0
if length(EbN0) > 1
    for j = 1:2
        for k in eachindex(EbN0)
            for i in 1:L
                str = curves[i,1]
                x = readdlm(directory*FB[j]*"ER_"*str*".txt",'\t',Float64,'\n')
                x = x[:,k]
                if j == 1
                    FER_EbN0[k,i] = log10.(x[iter])
                else
                    BER_EbN0[k,i] = log10.(x[iter])
                end
            end
        end
    end

    p = plot(
        EbN0,FER_EbN0,
        title = "FER x EbN0 $Protocol (N = $N, R = $(R[1])/$(R[2]), Iter = $iter)",
        color=permutedims(curves[:,4]),
        label=permutedims(curves[:,1]),
        markershape=permutedims(curves[:,2]),
        line=permutedims(curves[:,3]),
        lw=4,
        markersize=8,
        markerstrokewidth = 1,
        guidefontsize=10,
        legend_position = :bottomleft,
        tickfontsize=16,
        legend_font_pointsize = 20,
        size = (1400,1000)
    )
    display(p)
end
