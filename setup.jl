################################################################################
# Allan Eduardo Feitosa
# 26 mai 2025
# Setup of the LDPC simulation

############################## ) INCLUDED FILES ###############################

include("NR LDPC/NR_LDPC_parameters.jl")
include("IEEE80216e.jl")
include("PEG.jl")
include("LU_encoding.jl")
include("find_girth.jl")
include("make_VN_CN_lists.jl")
include("print_simulation_details.jl")
include("prepare_simulation.jl")

############################# PARITY-CHECK MATRIX #############################
if PROTOCOL == "NR5G"
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
    LL = nothing
    UU = nothing
else
    NN = GG
    AA = round(Int,GG*RR)
    KK, g_CRC = get_CRC_poly(AA)
    _, G_CRC = get_CRC_poly(AA)  
    if PROTOCOL == "PEG"     
        LIFTSIZE = 0
        E_H = nothing
        MM = NN - KK
        # Generate Parity-Check Matrix by the PEG algorithm
        H_PEG, GIRTH = PEG(LAMBDA,RO,MM,NN)
        HH, LL, UU = remake_H(H_PEG,0)
    elseif PROTOCOL == "WiMAX"
        # N takes values in {576,672,768,864,960,1056,1152,1248,1344,1440,1536,
        # 1632,1728,1824,1920,2016,2112,2208,2304}.    
        # R takes values in {"1/2","2/3A","2/3B","3/4A","3/4B","5/6"}.
        HH, LIFTSIZE, E_H = IEEE80216e(NN,RR,"A")
        MM,_ = size(HH)
        KK = NN - MM
        AA = KK - 16 # since max(AA) = 2304*5/6 ≤ 3824, g_cRC = {CRC16}      
        GIRTH = find_girth(HH,100000)
        LL = nothing
        UU = nothing
    end
end

# list of checks and variables nodes
NC = make_cn2vn_list(HH)
NV = make_vn2cn_list(HH)

STR = 
"""############################### LDPC parameters ################################
LDPC Protocol: """
if PROTOCOL == "NR5G"
    STR *= "NR-LDPC (5G), Zc = $(NR_LDPC_DATA.Zc)"
elseif PROTOCOL == "PEG"
    STR *= "PEG"
elseif PROTOCOL == "WiMAX"
    STR *= "IEEE80216e (WiMAX)"
end
STR *= "\nParity Check Matrix: $MM x $NN"
println(STR)
if SAVE
    println(FILE,STR)
end

display(sparse(HH))

STR = """

Graph girth = $GIRTH
Effective rate = $(round(RR,digits=3))
"""
println(STR)
if SAVE
    println(FILE,STR)
end

# Number of Threads
if !TEST
    NTHREADS = Threads.nthreads()
else
    NTHREADS = 1
end

####################### GENERATE NOISE AND SAMPLE SEEDS ########################

RGN_NOISE_SEEDS = zeros(Int,NTHREADS)
RGN_MESSAGE_SEEDS = zeros(Int,NTHREADS)
for i in eachindex(RGN_NOISE_SEEDS)
    RGN_NOISE_SEEDS[i] = SEED_NOISE + i - 1
    RGN_MESSAGE_SEEDS[i] = SEED_MESSA + 1 - 1
end

################################## SIMULATION ##################################
if TEST
    LRM = Dict()
    LQM = Dict()
    for i in eachindex(ACTIVE)
        if ACTIVE[i]
            if PROF
                global PRIN = false
                if TRIALS_TEST < 1000
                    global TRIALS_TEST = 1000
                end
                if MAXITER_TEST < 20
                    global MAXITER_TEST = 20
                end
                @profview _,_ = prepare_simulation(
                                                EbN0_TEST,
                                                MODES[i],
                                                TRIALS_TEST,
                                                MAXITER_TEST,
                                                BPTYPES[i],
                                                DECAY_TEST)
            else
            LRM[MODES[i]], LQM[MODES[i]]  = prepare_simulation(
                                                EbN0_TEST,
                                                MODES[i],
                                                TRIALS_TEST,
                                                MAXITER_TEST,
                                                BPTYPES[i],
                                                DECAY_TEST)
            end
        end
    end
else
    if STOP
        FER_LABELS = Vector{String}()
        FERMAX = Vector{Vector{<:AbstractFloat}}()
    end
    FER = Dict()
    BER = Dict()
    for i in eachindex(ACTIVE)
        if ACTIVE[i]
            for decay in DECAYS[i]
                if decay != 0.0
                    mode = MODES[i]*" $decay"
                else
                    mode = MODES[i]
                end
                FER[mode], BER[mode] = prepare_simulation(
                                        EbN0,
                                        MODES[i],
                                        TRIALS,
                                        MAXITERS[i],
                                        BPTYPES[i],
                                        decay)
                if STOP
                    push!(FER_LABELS,mode*" ($(BPTYPES[i]))")
                    push!(FERMAX,FER[mode][MAXITERS[i],:])
                end
            end
        end
    end
end  