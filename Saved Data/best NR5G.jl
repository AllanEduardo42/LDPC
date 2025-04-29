using DelimitedFiles
using Plots

decays = [0.7, 0.8, 0.9, 1.0, 0.2]

plots = [true,
        true,
        false,
        false,
        false,
        true,
        false,
        false,
        true,
        false,
        true,
        true,
        true
        ]

lines = [:solid,
         :solid,
         :dash,
         :dash,
         :dash,
         :dash,
         :dash,
         :dash,
         :dash,
         :dash,
         :dot,
         :dot,
         :dashdot
         ]

maxiter = 10

index = [3,
         3,
         3,
         2,
         3,
         2,
         3,
         2,
         3,
         1,
         2,
         5,
         3
         ]

plotlyjs()
directory1 = "./Saved Data/"
directory2 = ["RBP NR5G",
              "RBP NR5G relative",
              "List-RBP NR5G 128 16",
              "List-RBP NR5G 128 16 relative",
              "List-RBP NR5G 128 8",
              "List-RBP NR5G 128 8 relative (all decays)",
              "List-RBP NR5G 128 4",
              "List-RBP NR5G 128 4 relative",
              "List-RBP NR5G 128 2",
              "List-RBP NR5G 128 2 relative",
              "VN-RBP NR5G (all decays)",
              "VN-RBP NR5G relative (all decays)",
              "NS NR5G (all decays)"
              ]
file = ["FER_RBP",
        "FER_RBP",
        "FER_List-RBP",
        "FER_List-RBP",
        "FER_List-RBP",
        "FER_List-RBP",
        "FER_List-RBP",
        "FER_List-RBP",
        "FER_List-RBP",
        "FER_List-RBP",
        "FER_VN-RBP",
        "FER_VN-RBP",
        "FER_NS-RBP"
        ]

lim = log10(1/maximum(512000))

labels = ["RBP : ", 
          "RBP rel : ",
          "List (128,16) : ",
          "List rel (128,16) : ",
          "List (128,8) : ",
          "List rel (128,8) : ",
          "List (128,4) : ",
          "List rel (128,4) : ",
          "List (128,2) : ",
          "List rel (128,2) : ",
          "VN-RBP : ",
          "VN-RBP rel : ",
          "NS : "
          ]


title = "NR5G FER (SNR = 1.8 dB)"
p = plot()
for i in eachindex(plots)
    if plots[i]
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
            ylims=(lim,0),
            # xlims=(1,maxiter),
            # color=[1 2 3 4],
            legend_title = "Algorithm : Decay factor",
            legend_title_font_pointsize = 9,
            # legend = :outertopright,
            size = (900,600)
        )
    end
end
display(p)
# Plots.pdf(p,"rbp_opt_2.pdf")
