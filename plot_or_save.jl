################################################################################
# Allan Eduardo Feitosa
# 26 mai 2025
# Save or plot the results of the LDPC simulation

################################## SAVE DATA ###################################
if SAVE        
    if length(EbN0) > 1
        # FER
        open("./Saved Data/"*NOW*"/FERMAX.txt","w") do io
                        writedlm(io,FERMAX)
        end
        
        # BER
        open("./Saved Data/"*NOW*"/BERMAX.txt","w") do io
                        writedlm(io,BERMAX)
        end
    end
    if !STOP
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
    LIMFER = 1 ./(TRIALS)
    LIMBER = 1 ./((TRIALS)*NN)

    if !STOP  
        FB = ["F","B"]        
        for j=1:2
            for k in eachindex(EbN0)
                p = plot()
                title = FB[j]*"ER $PROTOCOL (EbN0 = $(EbN0[k])dB, R = $(round(RR,digits=2)), G = $GG)"
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
                                y = log10.(FER[mode][:,k])
                                lim = log10(LIMFER[k])
                            else
                                y = log10.(BER[mode][:,k])
                                lim = log10(LIMBER[k])
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
                                ylims=(lim,0),
                                size = (700,500)
                            )
                        end
                    end
                end
            display(p)
            end
        end
        # PLOT = plot!(1:MAXITERS[i],
        #       ones(Int,MAXITER)*LIMFER,
        #       label = "minimum",
        #       lw=2,
        #       color=:black,
        #       ls=:dash)
    end

    LABELS = permutedims(LABELS)
    
    if length(EbN0) > 1
        # FER x EbN0
        PLOT = plot(
            EbN0,FERMAX,
            xlabel="EbN0 (dB)",
            label=LABELS,
            lw=2,
            title="FER x EbN0",
            ylims=(log10(LIMFER[end]),0)
        )
        display(PLOT)
        # BER x EbN0
        PLOT = plot(
            EbN0,BERMAX,
            xlabel="EbN0 (dB)",
            label=LABELS,
            lw=2,
            title="BER x EbN0",
            ylims=(log10(LIMBER[end]),0)
        )
        display(PLOT)
    end
end
println("The End!")
if SAVE
    println(FILE, "The End!")
    close(FILE)
end