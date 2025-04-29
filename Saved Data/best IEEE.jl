using DelimitedFiles
using Plots

decays = [0.7, 0.8, 0.9, 1.0]

plots = [true,
        true,
        true,
        false,
        false,
        false,
        false,
        false,
        false,
        true,
        true,
        true,
        true,
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

index = [2,
         3,
         4,
         2,
         3,
         3,
         3,
         2,
         2,
         2,
         2,
         4,
         3
         ]

plotlyjs()
directory1 = "./Saved Data/"
directory2 = ["RBP IEEE (all decays)",
              "RBP IEEE relative (all decays)",
              "List-RBP IEEE 128 2 (all decays)",
              "List-RBP IEEE 128 2 (all decays) relative",
              "List-RBP IEEE 128 4 (all decays)",
              "List-RBP IEEE 128 4 (all decays) relative",
              "List-RBP IEEE 128 8 (all decays)",
              "List-RBP IEEE 128 8 (all decays) relative",
              "List-RBP IEEE 128 16 (all decays)",
              "List-RBP IEEE 128 16 (all decays) relative",
              "VN-RBP IEEE (all decays)",
              "VN-RBP IEEE (all decays) relative",
              "NS IEEE"
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
          "RBP rel: ", 
          "List RBP (128,2) : ",
          "List RBP (128,2) rel :",
          "List RBP (128,4) : ",
          "List RBP (128,4) rel :",
          "List RBP (128,8) : ",
          "List RBP (128,8) rel : ",
          "List RBP (128,16) : ",
          "List RBP (128,16) rel : ",
          "VN-RBP : ",
          "VN-RBP rel : ",
          "NS : "
          ]


title = "IEEE FER (SNR = 1.8 dB)"
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
            # xlims=(1,10),
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
