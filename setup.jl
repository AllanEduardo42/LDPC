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
include("get_CRC_poly.jl")
include("GF2_poly.jl")
include("find_girth.jl")
include("make_VN_CN_lists.jl")
include("print_simulation_details.jl")
include("print_algorithm_details.jl")
include("prepare_simulation.jl")
include("simcore.jl")

if PROTOCOL == "5GNR"
    include("5G NR LDPC/NR_LDPC_parameters.jl")
    include("5G NR LDPC/NR_LDPC_auxiliary_functions.jl")
    include("5G NR LDPC/NR_LDPC_make_parity_check_matrix.jl")
elseif PROTOCOL == "WiMAX"
    include("IEEE80216e.jl")
elseif PROTOCOL == "PEG"
    include("PEG.jl")
    include("LU_encoding.jl")
else
    throw(error(lazy"Invalid coding protocol (valid options: 5GNR, WiMAX and PEG)"))
end

# simcore functions
include("Simulation core functions/auxiliary_functions.jl")
include("Simulation core functions/encode_LDPC.jl")
include("Simulation core functions/calc_prior_LLRs.jl")
include("Simulation core functions/tanh_V2C.jl")
include("Simulation core functions/calc_syndrome.jl")
include("Simulation core functions/calc_C2V.jl")
include("Simulation core functions/calc_post_LLR.jl")

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
if PROTOCOL == "5GNR"
    H1 = [false false]
    RV = 0
    AA, KK, RR, G_CRC, LIFTSIZE, NR_LDPC_DATA = NR_LDPC_parameters(GG,RR,RV,false)
    HH, E_H = NR_LDPC_make_parity_check_matrix(LIFTSIZE,
                                               NR_LDPC_DATA.iLS,
                                               NR_LDPC_DATA.bg,
                                               NR_LDPC_DATA.P,
                                               NR_LDPC_DATA.K_prime,
                                               NR_LDPC_DATA.K,
                                               NR_LDPC_DATA.P_Zc)
    MM, NN = size(HH)
    GIRTH = find_girth(HH,100000)
    LL = H1
    UU = H1
else
    NN = GG
    AA = round(Int,GG*RR)
    KK, G_CRC = get_CRC_poly(AA)  
    if PROTOCOL == "PEG"     
        LIFTSIZE = 0
        E_H = [0 0]
        MM = NN - KK
        # Generate Parity-Check Matrix by the PEG algorithm
        H_PEG, GIRTH = PEG(LAMBDA,RO,MM,NN)
        HH, LL, UU = remake_H(H_PEG,0)
        H1 = HH[:,1:KK]
    elseif PROTOCOL == "WiMAX"
        # N takes values in {576,672,768,864,960,1056,1152,1248,1344,1440,1536,
        # 1632,1728,1824,1920,2016,2112,2208,2304}.    
        # R takes values in {"1/2","2/3A","2/3B","3/4A","3/4B","5/6"}.
        HH, LIFTSIZE, E_H = IEEE80216e(NN,RR,"A")
        MM,_ = size(HH)
        KK = NN - MM
        AA = KK - 16 # since max(AA) = 2304*5/6 ≤ 3824, g_cRC = {CRC16}      
        GIRTH = find_girth(HH,100000)
        LL = [false false]
        UU = [false false]
        H1 = [false false]
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

    print_simulation_details(ERRORS_TEST,MAXITER_TEST,EbN0_TEST)

    LRM = Dict()
    LQM = Dict()
    COUNT = 0
    for i in eachindex(ACTIVE)
        if ACTIVE[i]
            global COUNT += 1
            if DECAYS[i][1] == 0.0
                global DECAY_TEST = 0.0
            end
            for bptype in BPTYPES[i]
                if PROF
                    global PRIN = false
                    # evaluate expr for profile view
                    prog = "@profview (prepare_simulation([EbN0_TEST],
                            ALGORITHMS[i],
                            TRIALS_TEST,
                            MAXITER_TEST,
                            bptype,
                            DECAY_TEST))"
                    expr = Meta.parse(prog)
                    eval(expr)
                else
                LRM[ALGORITHMS[i]], LQM[ALGORITHMS[i]] = prepare_simulation(
                                                    COUNT,
                                                    [EbN0_TEST],
                                                    ALGORITHMS[i],
                                                    ERRORS_TEST,
                                                    MAXITER_TEST,
                                                    bptype,
                                                    DECAY_TEST)
                end
            end
        end
    end
else

    print_simulation_details(ERRORS,MAXITER,EbN0)

    FER = Dict()
    BER = Dict()
    COUNT = 0
    for i in eachindex(ACTIVE)
        if ACTIVE[i]
            global COUNT += 1
            algo = ALGORITHMS[i]
            # if algo == "List-RBP"
            #     algo *= " ($(LISTSIZES[1]),$(LISTSIZES[2]))"
            # end
            for bptype in BPTYPES[i]
                if bptype != "TANH"
                    algo *= " ($bptype)"
                end
                for decay in DECAYS[i]
                    if length(DECAYS[i]) > 1 && decay != 0.0
                        algo *= " $decay"
                    end
                    fer, ber = prepare_simulation(
                                            COUNT,
                                            EbN0,
                                            ALGORITHMS[i],
                                            ERRORS,
                                            MAXITERS[i],
                                            bptype,
                                            decay)
                    if SAVE
                        open("./Saved Data/"*NOW*"/FER_"*algo*".txt","w") do io
                                writedlm(io,fer)
                        end
                        open("./Saved Data/"*NOW*"/BER_"*algo*".txt","w") do io
                                writedlm(io,ber)
                        end
                    end
                    FER[algo] = fer
                    BER[algo] = ber
                end
            end
        end
    end
end  