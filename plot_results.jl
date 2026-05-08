################################################################################
# Allan Eduardo Feitosa
# 28 ago 2025
# Plot the results of the LDPC simulation

MARKERS = [:none,:none,:none,:dtriangle,:circle,:rect,:utriangle,:diamond,:cross,
           :star5,:hexagon,:none]

plotlyjs()
# gr()

# LIMFER = 1 ./(TRIALS)
# LIMBER = 1 ./((TRIALS)*NN)

LABELS = Vector{String}()
FERMAX = zeros(length(EbN0),NUM_MODES)
BERMAX = zeros(length(EbN0),NUM_MODES)

 
FB = ["F","B"]        
for j=1:2          
    for k in eachindex(EbN0)  
        p = plot()      
        title = FB[j]*"ER $PROTOCOL (Code length = $CODE_LENGTH, Rate = $(RATE[1])/$(RATE[2]), Eb/N0 = $(EbN0[k])dB)"      
        for i in eachindex(ACTIVE)
            if ACTIVE[i]
                algo = ALGORITHMS[i]
                for decay in DECAYS[i]
                    if length(DECAYS[i]) > 1 && decay != 0.0
                        algo *= " $decay"
                    end
                    if j == 1 && k == 1
                        push!(LABELS,algo)
                    end
                    if FB[j] == "F"
                        y = log10.(FER[algo][:,k])
                    else
                        y = log10.(BER[algo][:,k])
                    end
                    x = 1:MAXITER
                    p = plot!(
                        x,
                        y,
                        xlabel="Iteration",
                        label=algo,
                        lw=2,
                        title=title,
                        size = (800,500)
                    )
                end
            end
        end
        display(p)
    end
end

LABELS = permutedims(LABELS)

TITLE = "ER $PROTOCOL (Code length = $CODE_LENGTH, Rate = $(RATE[1])/$(RATE[2]), Iter = $MAXITER)"

if length(EbN0) > 1

    # FER x EbN0
    for j=1:2
        p = plot()
        count = 0
        for i in eachindex(ACTIVE)
            if ACTIVE[i]
                count += 1
                if j == 1
                    y = FERMAX[:,i]
                else
                    y = BERMAX[:,i]
                end
                p = plot!(
                    EbN0,y,
                    label=LABELS[count],
                    lw=2,
                    title=FB[j]*TITLE,
                    size = (700,500)
                )
            end
        end
        display(p)
    end
end

println("The End!")