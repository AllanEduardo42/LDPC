################################################################################
# Allan Eduardo Feitosa
# 26 mai 2025
# Setup of the LDPC simulation

if ACTIVE_ALL
    for i in eachindex(ACTIVE)
        if ACTIVE[i]
            ACTIVE[i] = false
        else
            ACTIVE[i] = true
        end
    end
end

################################ INCLUDED FILES ################################

include("GF2_functions.jl")
include("auxiliary setup functions.jl")
include("print_simulation_details.jl")
include("print_algorithm_details.jl")
include("prepare_simulation.jl")
include("simcore.jl")

if PROTOCOL == "5GNR"
    include("5G NR LDPC/NR_LDPC_functions.jl")
    include("5G NR LDPC/NR_LDPC_parameters.jl")
    include("5G NR LDPC/NR_LDPC_auxiliary_functions.jl")
    include("5G NR LDPC/NR_LDPC_make_parity_check_matrix.jl")
elseif PROTOCOL == "WiMAX"
    include("WiMAX.jl")
elseif PROTOCOL == "PEG"
    include("PEG.jl")
    include("LU_encoding.jl")
else
    throw(error(lazy"Invalid coding protocol (valid options: 5GNR, WiMAX and PEG)"))
end

# algorithms
for i in eachindex(ACTIVE)
    if ACTIVE[i]
        if ALGORITHMS[i] == "C-RBP" || ALGORITHMS[i] == "C&DR-RBP"
            include("Algorithms/C&R-RBP.jl")
        else
            include("Algorithms/$(ALGORITHMS[i]).jl")
        end
    end
end

############################# PARITY-CHECK MATRIX #############################

RR = RATE[1]/RATE[2]

if PROTOCOL == "5GNR"
    RV = 0
    AA, KK, RR, G_CRC, LIFTSIZE, NR_LDPC_DATA = NR_LDPC_parameters(CODE_LENGTH,RR,RV,false)
    HH, E_H = NR_LDPC_make_parity_check_matrix(LIFTSIZE,
                                               NR_LDPC_DATA.iLS,
                                               NR_LDPC_DATA.bg,
                                               NR_LDPC_DATA.P,
                                               NR_LDPC_DATA.K_prime,
                                               NR_LDPC_DATA.K,
                                               NR_LDPC_DATA.P_Zc)
    MM, NN = size(HH)
    GIRTH = find_girth(HH,100000)
    PP = nothing
else
    NN = CODE_LENGTH
    G_CRC = nothing
    if PROTOCOL == "PEG"     
        KK = round(Int,CODE_LENGTH*RR)
        LIFTSIZE = 0
        E_H = [0 0]
        MM = NN - KK
        # Generate Parity-Check Matrix by the PEG algorithm
        H_PEG, GIRTH = PEG(LAMBDA,RO,MM,NN)
        HH, LL, UU = LU_encoding(H_PEG,0)
        PP = generate_parity_matrix(HH,LL,UU)        
    elseif PROTOCOL == "WiMAX"
        # N takes values in {576,672,768,864,960,1056,1152,1248,1344,1440,1536,
        # 1632,1728,1824,1920,2016,2112,2208,2304}.    
        # R takes values in {"1/2","2/3A","2/3B","3/4A","3/4B","5/6"}.
        HH, LIFTSIZE, E_H = WiMAX(NN,RR,"A")
        MM,_ = size(HH)
        KK = NN - MM   
        GIRTH = find_girth(HH,100000)
        PP = nothing
    end
end

# list of checks and variables nodes
NC = make_cn2vn_list(HH)
NV = make_vn2cn_list(HH)

# Number of Threads
if !TEST
    NTHREADS = Threads.nthreads()
else
    NTHREADS = 1
end

####################### GENERATE NOISE AND SAMPLE SEEDS ########################

RGN_SEEDS = zeros(Int,NTHREADS)
for i in 1:NTHREADS
    RGN_SEEDS[i] = SEED + i - 1
end

############################## PREPARE SIMULATION ##############################
if TEST

    print_simulation_details(MAX_FRAME_ERRORS_TEST,MAXITER_TEST,EbN0_TEST)

    LRM = Dict()
    LQM = Dict()
    COUNT = 0
    for i in eachindex(ACTIVE)
        if ACTIVE[i]
            global COUNT += 1
            if !is_RD_algorithm(ALGORITHMS[i])
                global DECAY_TEST = 0.0
            end
            if PROF
                global PRIN = false
                # evaluate expr for profile view
                global ALGO_PROF = ALGORITHMS[i]
                prog = "@profview (prepare_simulation(
                        COUNT,
                        [EbN0_TEST],
                        ALGO_PROF,
                        MAX_FRAME_ERRORS_TEST,
                        MAXITER_TEST,
                        BPTYPE,
                        DECAY_TEST))"
                expr = Meta.parse(prog)
                eval(expr)
            else
            LRM[ALGORITHMS[i]], LQM[ALGORITHMS[i]] = prepare_simulation(
                                                COUNT,
                                                [EbN0_TEST],
                                                ALGORITHMS[i],
                                                MAX_FRAME_ERRORS_TEST,
                                                MAXITER_TEST,
                                                DECAY_TEST)
            end
        end
    end
else

    print_simulation_details(MAX_FRAME_ERRORS,MAXITER,EbN0)

    FER = Dict()
    BER = Dict()
    LABELS = []
    COUNT = 0
    for i in eachindex(ACTIVE)
        if ACTIVE[i]
            if is_RD_algorithm(ALGORITHMS[i])
                decays = DECAYS
            else
                decays = [0.0]
            end
            for decay in decays
                global COUNT += 1
                algo = ALGORITHMS[i]
                if length(decays) > 1
                    algo *= " $decay"                  
                end
                fer, ber = prepare_simulation(
                                        COUNT,
                                        EbN0,
                                        ALGORITHMS[i],
                                        MAX_FRAME_ERRORS,
                                        MAXITER,
                                        decay)
                if SAVE
                    open(DIRECTORY*"/FER_"*algo*".txt","w") do io
                        writedlm(io,fer)
                    end
                    open(DIRECTORY*"/BER_"*algo*".txt","w") do io
                        writedlm(io,ber)
                    end
                end
                FER[algo] = fer
                BER[algo] = ber
                push!(LABELS,algo)
            end
        end
    end
end  