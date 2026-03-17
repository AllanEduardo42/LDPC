using DelimitedFiles
using Plots

N = 512
R = [1,2]
EbNo = 1.5
protocol = "NR5G"
decays = collect(0.1:0.1:1.0)
FB = ["F","B"]
maxiter = 50
# index = [1,2,2,1,2]

best = 7

plotlyjs()
directory = "./Saved Data/VN-RBP $protocol $N $(R[1])|$(R[2]) $EbNo all decays/"
lim = log10(1/maximum(1000*2^10))+1

for j=1:2
    title = FB[j]*"ER $protocol (N = $N, R = $(R[1])/$(R[2]), Eb/N0 = $(EbNo)dB)"
    p = plot()
    sdecays = string.(decays)
    for i in eachindex(sdecays)
        x = readdlm(directory*FB[j]*"ER_VN-RBP "*sdecays[i]*".txt",'\t',Float64,'\n')
        x = x[1:maxiter,:]
        labels = Vector{String}()
        push!(labels,"VN-RBP "*sdecays[i])
        labels = permutedims(labels)
        if i == best
            line = :solid
        else
            line = :dot
        end
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
            color=i,
            legend_title = "(decay, EbNo)",
            legend_title_font_pointsize = 8,
            legend = :outertopright,
            size = (900,600)
        )

    end
    display(p)
    # Plots.pdf(p,"rbp.pdf")
    global lim -= 1
end

