using DelimitedFiles
using Plots

decay = 0.9
lines = [:solid,:solid,:dash,:dot,:dot,:solid]
maxiter = 30

index = [3,3,3,2,3,4,3,3,3,3]

plotlyjs()
directory1 = "./Saved Data/"
directory2 = ["RBP decay test","RBP decay test relative",
              "RBP decay test new 1", "RBP decay test new 2",
              "RBP decay test new 3", "RBP decay test new 4"]

file = "FER_RBP"
lim = log10(1/maximum(1000*2^10))-0.5

labels = ["RBP", "RBP relative 0",
         "RBP relative 1", "RBP relative 2",
         "RBP relative 3: ", "RBP relative 4"]


title = "FER (SNR = 2.0 dB, decay = $decay) "
p = plot()
for i in eachindex(directory2)
    x = readdlm(directory1*directory2[i]*"/"*file*" "*string(decay)*".txt",'\t',Float64,'\n')
    display(x)
    p = plot!(
        1:maxiter,
        x[:,4],
        xlabel="Iteration",
        label=labels[i],
        lw=2,
        ls=lines[i],
        title=title,
        ylims=(lim,0),
        xlims=(0,30),
        # color=[1 2 3 4],
        legend_title = "Algorithm : Decay factor",
        legend_title_font_pointsize = 9,
        # legend = :outertopright,
        size = (900,600)
    )
end
display(p)
Plots.pdf(p,"rbp_opt_2.pdf")
