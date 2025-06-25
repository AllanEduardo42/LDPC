using DelimitedFiles
using Plots

function save_pdf(p, filename)
    width, height = Int.(p.attr[:size])
    Plots.prepare_output(p)
    PlotlyJS.savefig(Plots.plotlyjs_syncplot(p), filename*".pdf", width=width, height=height)
end

plotlyjs()

FB = ["F","B"]
markers = [:none, :none, :circle, :rect, :utriangle, :diamond, :cross]
modes = ["Flooding","LBP","RBP","RBP relative", "List-RBP","NW-RBP","VN-RBP"]
directory = "./Saved Data/Testes $protocol $N $(R[1])|$(R[2]) $EbNo/"
liminf = log10(1/iters)+1
limsup = 0.1

for j=1:2
    title = FB[j]*"ER $protocol (N = $N, R = $(R[1])/$(R[2]), Eb/N0 = $(EbNo)dB)"
    p = plot()
    for k in eachindex(modes)
        if modes[k] == "Flooding" || modes[k] == "LBP"
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
            elseif modes[k] == "LBP"
                push!(labels,"LBP")
            elseif modes[k] == "List-RBP"
                push!(labels,modes[k]*" (16,2), d = "*sdecays[i])
            else
                push!(labels,modes[k]*", d = "*sdecays[i])
            end
            line = :solid
            if i == best[k]
                line = :solid
                marker = markers[k]
            else
                line = :dot
                marker = :none
            end
            labels = permutedims(labels)
            if !onlybest || i == best[k]
                p = plot!(
                    1:maxiter,
                    log10.(x),
                    xlabel="Iteration",
                    label=labels,
                    lw=2,
                    ls=line,
                    title=title,
                    ylims=(liminf,limsup),
                    # xlim=(1,maxiter),
                    color=k,
                    legend_title = "(algorithm, best decay factor)",
                    legend_title_font_pointsize = 10,
                    legend_font_pointsize = 10,
                    # legend = :outertopright,
                    size = 1.6 .*(600,400),
                    markershape=marker
                )
            end
        end
    end
    display(p)
    save_pdf(p,directory*"/$(FB[j])ER")
    global liminf -= 1.2
    global limsup -= 1
end



