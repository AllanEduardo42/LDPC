################################################################################
# Allan Eduardo Feitosa
# 26 mai 2025
# Save or plot the results of the LDPC simulation

################################## SAVE DATA ###################################
if SAVE        
    if STOP
        aux = []
        for i in eachindex(FERMAX)
            push!(aux,(FER_LABELS[i],FERMAX[i]))
        end
        FERS = Dict(aux)
        CSV.write("./Saved Data/"*NOW*"/FERMAX.csv", DataFrame(FERS), header=true)
    else
        for i in eachindex(ACTIVE)
            if ACTIVE[i]
                for decay in DECAYS[i]
                    if decay != 0.0
                        mode = MODES[i]*" $decay"
                    else
                        mode = MODES[i]
                    end
                    open("./Saved Data/"*NOW*"/FER_"*mode*".txt","w") do io
                        writedlm(io,FER[mode])
                    end
                    open("./Saved Data/"*NOW*"/BER_"*mode*".txt","w") do io
                        writedlm(io,BER[mode])
                    end
                end
            end
        end
    end
else
################################### PLOTTING ###################################
    plotlyjs()
    # gr()
    LIMFER = 1/maximum(TRIALS)
    LIMBER = 1/(maximum(TRIALS)*NN)

    # FER x EbN0
    if STOP && length(EbN0) > 1
        FER_LABELS = permutedims(FER_LABELS)
        PLOT = plot(
            EbN0,FERMAX,
            xlabel="EbN0 (dB)",
            label=FER_LABELS,
            lw=2,
            title="FER MAXITER",
            ylims=(LIMFER,0)
        )
        display(PLOT)
    end

    if !STOP  
        FB = ["F","B"]  
        # FER x Iterations
        for j=1:2
            p = plot()
            title = FB[j]*"ER $PROTOCOL (R = $(round(RR,digits=2)))"
            for i in eachindex(ACTIVE)
                if ACTIVE[i]
                    for decay in DECAYS[i]
                        if decay != 0.0
                            mode = MODES[i]*" $decay"
                        else
                            mode = MODES[i]
                        end
                        labels = Vector{String}()
                        for enr in EbN0
                            push!(labels,"$mode (Eb/N0 = $(enr)dB)")
                        end
                        labels = permutedims(labels)
                        if FB[j] == "F"
                            y = log10.(FER[mode])
                            lim = log10(LIMFER)
                        else
                            y = log10.(BER[mode])
                            lim = log10(LIMBER)
                        end
                        x = 1:MAXITERS[i]
                        p = plot!(
                            x,
                            y,
                            # yscale=:log10,
                            xlabel="Iteration",
                            # minorgrid=true,
                            label=labels,
                            lw=2,
                            title=title,
                            ylims=(lim,0)
                        )
                    end
                end
            end
            display(p)
        end
        # PLOT = plot!(1:MAXITERS[i],
        #       ones(Int,MAXITER)*LIMFER,
        #       label = "minimum",
        #       lw=2,
        #       color=:black,
        #       ls=:dash)
    end
end
println("The End!")
if SAVE
    println(FILE, "The End!")
    close(FILE)
end