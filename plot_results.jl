################################################################################
# Allan Eduardo Feitosa
# 28 ago 2025
# Plot the results of the LDPC simulation

plotlyjs()

FER_EbN0 = zeros(length(EbN0),COUNT)
BER_EbN0 = zeros(length(EbN0),COUNT)
 
FB = ["F","B"]        
       
for ebn0 in eachindex(EbN0)  
    for j = 1:2   
        p = plot()      
        title = FB[j]*"ER $PROTOCOL (Code length = $CODE_LENGTH, Rate = $(RATE[1])/$(RATE[2]), Eb/N0 = $(EbN0[ebn0])dB)"      
        for i in 1:COUNT
            algo = LABELS[i]
            if j == 1
                y = log10.(FER[algo][:,ebn0])
                FER_EbN0[ebn0,i] = y[ITER]
            else
                y = log10.(BER[algo][:,ebn0])
                BER_EbN0[ebn0,i] = y[ITER]
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
        display(p)
    end
end

LABELS = permutedims(LABELS)

TITLE = "ER $PROTOCOL (Code length = $CODE_LENGTH, Rate = $(RATE[1])/$(RATE[2]), Iter = $ITER)"

if length(EbN0) > 1

    # FER x EbN0
    for j=1:2
        if j == 1
            y = FER_EbN0
        else
            y = BER_EbN0
        end
        p = plot(
            EbN0,y,
            label=LABELS,
            lw=2,
            title=FB[j]*TITLE,
            size = (700,500)
        )
        display(p)
    end
end

println("The End!")