################################################################################
# Allan Eduardo Feitosa
# 28 ago 2025
# Plot the results of the LDPC simulation

plotlyjs()
# gr()

LIMFER = 1 ./(TRIALS)
LIMBER = 1 ./((TRIALS)*NN)

LABELS = Vector{String}()
FERMAX = zeros(length(EbN0),NUM_MODES)
BERMAX = zeros(length(EbN0),NUM_MODES)

if !STOP  
    FB = ["F","B"]        
    for j=1:2          
        for k in eachindex(EbN0)  
            p = plot()      
            title = FB[j]*"ER $PROTOCOL (G = $GG, R = $(Rational(RR)), Eb/N0 = $(EbN0[k])dB)"      
            for i in eachindex(ACTIVE)
                if ACTIVE[i]
                    mode = MODES[i]
                    if mode == "List-RBP"
                        mode *= " ($(LISTSIZES[1]),$(LISTSIZES[2]))"
                    end
                    for decay in DECAYS[i]
                        mode2 = mode
                        if decay != 0.0
                            mode2 *= " $decay"
                        end
                        labels = Vector{String}()
                        push!(LABELS,mode2)
                        if FB[j] == "F"
                            y = log10.(FER[mode2][:,k])
                            lim = log10(LIMFER[k])
                        else
                            y = log10.(BER[mode2][:,k])
                            lim = log10(LIMBER[k])
                        end
                        x = 1:MAXITERS[i]
                        if j == 1
                            FERMAX[k,i] = y[MAXITER]
                        else
                            BERMAX[k,i] = y[MAXITER]
                        end
                        p = plot!(
                            x,
                            y,
                            xlabel="Iteration",
                            label=mode2,
                            lw=2,
                            title=title,
                            ylims=(lim,0),
                            size = (700,500)
                        )
                    end
                end
            end
            display(p)
        end
    end
end

LABELS = permutedims(LABELS)

TITLE = "ER $PROTOCOL (G = $GG, R = $(Rational(RR)), Iter = $MAXITER)"

if length(EbN0) > 1

    # FER x EbN0
    for j=1:2
        p = plot()
        for i in eachindex(ACTIVE)
            if ACTIVE[i]
                if j == 1
                    y = FERMAX[:,i]
                else
                    y = BERMAX[:,i]
                end
                p = plot!(
                    EbN0,y,
                    label=LABELS,
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
if SAVE
    println(FILE, "The End!")
    close(FILE)
end