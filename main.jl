############################# PARITY-CHECK MATRIX #############################

Msg = rand(Xoshiro(SEED_MESSA),Bool,A)

if LDPC == 1
    rv = 0
    H, Cword, Zc = NR_LDPC_encode(Msg,R,rv)
    M::Int, N::Int = size(H)
    girth = "?"
else
    Zc = 0
    L = round(Int,A/R)
    if LDPC == 2
        N::Int = L           
        M::Int = L - A
        # Vector of the variable node degrees
        D = rand(Xoshiro(SEED_GRAPH),densities,N-M)
        # Generate Parity-Check Matrix by the PEG algorithm
        H, girth = PEG(D,M,N)
    elseif LDPC == 3
        H = IEEE80216e(L,R)
        M::Int,N::Int = size(H)
        Msg = Msg[1:M]
        L = N
        girth = "?"
    else
        throw(ArgumentError(
                    lazy"CHECK = $CHECK, but must be 1, 2 or 3"
                ))
    end
    H1 = H[:,1:N-M]
    w = H1*Msg
    Cword = [Msg;w]
    display(iszero(H*Cword))
end

# Number of Threads
NTHREADS::Int = min(Threads.nthreads(),TRIALS)

# Aggregate mode string name and number of iterations
modes = [(FLOO,"Flooding",MAX),
        (_LBP,"LBP",MAX),
        (iLBP,"iLBP",MAX),
        (_RBP,"RBP",MAXRBP),
        (RRBP,"Random-RBP",MAXRBP),
        (LRBP,"Local-RBP",MAXRBP),
        (LIST,"List-RBP",MAXRBP)]

# RBP decay contants
decay = Dict(modes[4][2] => DECAYRBP,
            modes[5][2] => DECAYRRBP,
            modes[6][2] => DECAYLRBP,
            modes[7][2] => DECAYLIST)

####################### GENERATE NOISE AND SAMPLE SEEDS ########################

Rgn_noise_seeds = zeros(Int,NTHREADS)
rgn_samples_seeds = zeros(Int,NTHREADS)
for i in eachindex(Rgn_noise_seeds)
    Rgn_noise_seeds[i] = SEED_NOISE + i - 1
    rgn_samples_seeds[i] = SEED_SAMPL + i - 1
end

######################### PRIN INFORMATION ON SCREEN ##########################
println()
print("############################### LDPC parameters #######################")
println("#########")
println()
println("Parity Check Matrix: $M x $N")
println()
display(sparse(H))
println()
println("Graph girth = ", girth)
println()
println("Msg (L = $(length(Msg))):")
for i in eachindex(Msg)
    print(Int(Msg[i]))
    if i%80 == 0
        println()
    end
end
println()
println()
println("Cword (L = $(length(Cword))):")
for i in eachindex(Cword)
    print(Int(Cword[i]))
    if i%80 == 0
        println()
    end
end
println()
println()

#################### JULIA COMPILATION (FOR SPEED) AND TEST ####################
LR = Dict()
LQ = Dict()
p = (TRIALS ≤ 2) ? PRIN : false
for mode in modes
    if mode[1]
        LR[mode[2]] , LQ[mode[2]] = performance_sim(
            SNRTEST,
            mode[2],
            min(TRIALS,2),
            mode[3],;
            printtest=p)
    end
end                             
############################ PERFORMANCE SIMULATION ############################
if TRIALS > 2
    fer_labels = Vector{String}()
    fermax = Vector{Vector{<:AbstractFloat}}()
    FER = Dict()
    BER = Dict()
    for mode in modes
        if mode[1]
            @time FER[mode[2]], BER[mode[2]] = performance_sim(
                SNR,
                mode[2],
                TRIALS,
                mode[3];
                rgn_samples_seeds=rgn_samples_seeds)

            push!(fer_labels,mode[2])
            push!(fermax,FER[mode[2]][mode[3],:])
        end            
    end
end
################################### PLOTTING ###################################
if TRIALS > 2
    plotlyjs()
    lim = log10(1/TRIALS)

    # FER x SNR
    if length(SNR) > 1
        fer_labels = permutedims(fer_labels)    
        p = plot(
            SNR,fermax,
            xlabel="SNR (dB)",
            label=fer_labels,
            lw=2,
            title="FER (Graph girth = $girth)",
            ylims=(lim,0)
        )
        PLOT ? display(p) : nothing
        SAVE ? savefig(p, "FER.png") : nothing
    end

    # BER x Iterations
    if !STOP

        # FER x Iterations
        labels = Vector{String}()
        for snr in SNR
            push!(labels,"SNR (dB) = $snr")
        end
        labels = permutedims(labels)
        for mode in modes
            if mode[1]
                if mode[2] == "RBP" || mode[2] == "Random-RBP" ||
                mode[2] == "Local-RBP" || mode[2] == "List-RBP"
                    titlefer = "FER $(mode[2]) (decay factor = $(decay[mode[2]]))"
                else
                    titlefer = "FER $(mode[2])"
                end
                local p = plot(
                    1:mode[3],
                    FER[mode[2]],                
                    xlabel="Iteration",
                    label=labels,
                    lw=2,
                    title=titlefer,
                    ylims=(lim,0)
                )
                PLOT ? display(p) : nothing
                SAVE ? savefig(p,"FER_"*mode[2]*".png") : nothing
            end
        end

        for mode in modes
            if mode[1]
                if mode[2] == "RBP" || mode[2] == "Random-RBP" || 
                   mode[2] == "Local-RBP" || mode[2] == "List-RBP"
                    titleber = "BER $(mode[2]) (decay factor = $(decay[mode[2]]))"
                else
                    titleber = "BER $(mode[2])"
                end
                local p = plot(
                    1:mode[3],
                    BER[mode[2]],                
                    xlabel="Iteration",
                    label=labels,
                    lw=2,
                    title=titleber,
                    ylims=(lim-2,0)
                )
                PLOT ? display(p) : nothing
                SAVE ? savefig(p,"BER_"*mode[2]*".png") : nothing
            end
        end
    end
end
################################### SAVE DATA ##################################
if TRIALS > 2 && SAVE

    aux = []
    for i in eachindex(fermax)
        push!(aux,(fer_labels[i],fermax[i]))
    end
    FERS = Dict(aux)
    CSV.write("FERMAX.csv", DataFrame(FERS), header=true)

    # aux = []
    # padding = zeros(MAX-MAXRBP)
    # for mode in modes
    #     if mode[1]
    #         for i in eachindex(SNR)
    #             title = mode[2] * " (SNR=$(SNR))"
    #             push!(aux,(title,BER[mode[2]][:,i]))
    #         end
    #     end
    # end
    # BERS = Dict(aux)
    # CSV.write("BERS.csv", DataFrame(BERS), header=true)

end