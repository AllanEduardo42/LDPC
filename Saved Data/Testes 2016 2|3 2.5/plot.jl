using DelimitedFiles
using Plots

EbNo = [1.5]
decays = [0.8, 0.9, 1.0]
lines = [:dot,:dash,:solid]
colors = [1, 2, 3, 4, 5, 6]
FB = ["F","B"]
maxiter = 50

modes = ["Flooding","RBP","RBP relative", "List-RBP","NW-RBP","VN-RBP" ]

plotlyjs()
directory = "./Saved Data/Testes 2016 2|3 2.5/"
lim = log10(1/maximum(1000*2^10))+1

for j=1:2
    title = FB[j]*"ER NR5G (N = 2016, R = 2/3, Eb/N0 = 2.5)"
    p = plot()
    for k in eachindex(modes)
        if modes[k] == "Flooding"
            sdecays = [""]
            pad = ""
        else
            sdecays = string.(decays)
            pad = " "
        end
        for i in eachindex(sdecays)
            x = readdlm(directory*FB[j]*"ER_"*modes[k]*pad*sdecays[i]*".txt",'\t',Float64,'\n')
            x = x[1:maxiter,:]
            labels = Vector{String}()
            if modes[k] == "Flooding"
                push!(labels,"Flooding")
                line = lines[3]
            else
                push!(labels,modes[k]*" "*sdecays[i])
                line = lines[i]
            end
            labels = permutedims(labels)
            p = plot!(
                1:maxiter,
                log10.(x),
                xlabel="Iteration",
                label=labels,
                lw=2,
                ls=line,
                title=title,
                ylims=(lim,0),
                xlim=(1,maxiter),
                color=colors[k],
                legend_title = "(decay, EbNo)",
                legend_title_font_pointsize = 8,
                legend = :outertopright,
                size = (900,600)
            )
        end
    end
    display(p)
    # Plots.pdf(p,"rbp.pdf")
    global lim -= 1
end

