using DelimitedFiles
using Plots

decays = [0.7, 0.8, 0.9, 1.0]
lines = [:solid,
         :solid,
        #  :solid,
        #  :solid,
        #  :dash,
        #  :dash,
         :dot,
         :dot,
        #  :dashdot,
        #  :dashdot,
         :dash
         ]

maxiter = 30

index = [3,
         3,
        #  3,
        #  2,
        #  3,
        #  4,
         3,
         3,
        #  3,
        #  3,
         2
         ]

plotlyjs()
directory1 = "./Saved Data/"
directory2 = ["RBP decay test",
              "RBP decay test relative",
            #   "Local-RBP decay test",
            #   "Local-RBP decay test relative",
            #   "List-RBP decay test listsize2 = 1",
            #   "List-RBP decay test listsize2 = 1 relative",
              "List-RBP decay test listsize2 = 2",
              "List-RBP decay test listsize2 = 2 relative",
            #   "List-RBP decay test listsize2 = 4",
            #   "List-RBP decay test listsize2 = 4 relative",
              "VN-RBP test"]
file = ["FER_RBP",
        "FER_RBP",
        # "FER_Local-RBP",
        # "FER_Local-RBP",
        # "FER_List-RBP",
        # "FER_List-RBP",
        "FER_List-RBP",
        "FER_List-RBP",
        # "FER_List-RBP", 
        # "FER_List-RBP",
        "FER_VN-RBP"]

lim = log10(1/maximum(1000*2^10))

labels = ["RBP : ", 
          "RBP rel : ",
        #   "Local : ",
        #   "Local rel : ",
        #   "List (16,1) : ",
        #   "List rel (16,1) : ",
          "List (16,2) : ",
          "List rel (16,2) : ",
        #   "List (16,4) : ",
        #   "List rel (16,4) : ",
          "VN-RBP : "]


title = "FER (SNR = 2.0 dB)"
p = plot()
for i in eachindex(directory2)
    x = readdlm(directory1*directory2[i]*"/"*file[i]*" "*string(decays[index[i]])*".txt",'\t',Float64,'\n')
    x = x[1:maxiter,:]
    display(x)
    p = plot!(
        1:maxiter,
        x[:,4],
        xlabel="Iteration",
        label=labels[i]*" "*string(decays[index[i]]),
        lw=2,
        ls=lines[i],
        title=title,
        ylims=(lim,-2),
        xlims=(1,15),
        # color=[1 2 3 4],
        legend_title = "Algorithm : Decay factor",
        legend_title_font_pointsize = 9,
        # legend = :outertopright,
        size = (900,600)
    )
end
display(p)
# Plots.pdf(p,"rbp_opt_2.pdf")
