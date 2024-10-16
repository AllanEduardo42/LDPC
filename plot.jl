################################### PLOTTING ###################################
using Plots
using CSV, DataFrames

TRIALS::Int = 640320

SNR = collect(1:1:4)

FERS = CSV.read("FERS_RD.csv", DataFrame)

fer_labels = names(FERS)
fermax = Matrix(FERS)

if TRIALS > 1
    plotlyjs()
    lim = log10(1/TRIALS)
    # FER x SNR
    fer_labels = permutedims(fer_labels)    
    p = plot(
        SNR,fermax,
        xlabel="SNR (dB)",
        label=fer_labels,
        lw=2,
        title="FER (Graph girth = 8)",
        ylims=(lim,0)
    )
    display(p)
end