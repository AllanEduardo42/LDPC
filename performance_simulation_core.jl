################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Core routine to estimate the LPCD performance (FER BER x SNR)

include("auxiliary_functions.jl")
include("lookupTable.jl")
include("BP.jl")
include("calc_Lf.jl")

function
    performance_simulation_core(
        c::Vector{Bool},
        snr::Real,
        H::BitMatrix,
        mode::String,
        supermode::String,
        nreals::Integer,
        max::Integer,
        stop::Bool,
        rbpfactor::Union{AbstractFloat,Nothing},
        rgn_seed_noise::Integer,
        rng_seed_sample::Integer,
        test::Bool,
        t_test::Union{Vector{<:AbstractFloat},Nothing},
        printing::Bool
    )
    
################################## CONSTANTS ###################################
    
    # transform snr in standard deviations
    variance = 1 ./ (exp10.(snr/10))
    stdev = sqrt.(variance)
    # BPKS
    u = Float64.(2*c .- 1)
    # list of checks and variables nodes
    vn2cn  = make_vn2cn_list(H)
    cn2vn  = make_cn2vn_list(H)

    ############################# PREALLOCATIONS ################################

    # frame error rate
    FER = 0

    # bit error rate
    BER = zeros(max)
    ber = zeros(max)    

    # iteration in which SPA stopped
    # iters = Vector{Int}(undef, nreals)

    # estimate
    d = zeros(Bool,N)

    # syndrome
    syndrome = ones(Bool,M)

    # prior llr (if mode == "MKAY" just the prior probabilities)
    Lf = (mode != "MKAY") ? Vector{Float64}(undef,N) : Matrix{Float64}(undef,N,2)

    # noise
    noise = Vector{Float64}(undef,N)

    # received signal
    t = Vector{Float64}(undef,N)

    # bit-error
    bit_error = Vector{Bool}(undef,N) 

    # Lq -> matrix of vn2cn messages (N x M)
    # Lr -> matrix of cn2vn messages (M x N)
    # if mode == "MKAY" the are matrices for bit = 0 and bit = 1
    Lq, Lr = (mode != "MKAY") ? (H'*0.0,H*0.0) : (zeros(N,M,2),zeros(M,N,2))

    Ms = (supermode == "RBP") ? H*0.0 : nothing

    # Set variables Lrn and signs depending on the mode (also used for dispatch)
    if (mode == "TANH" && FAST) || mode == "LBP" || mode == "iLBP"
        Lrn = zeros(N) 
        signs = nothing
    elseif mode == "ALTN" || mode == "TABL"
        Lrn = zeros(N)
        signs = zeros(Bool,N)
    elseif mode == "MSUM" || supermode == "RBP"
        Lrn = nothing
        signs = zeros(Bool,N)
    else
        Lrn = nothing
        signs = nothing
    end

   # Set other variables that depend on the mode

    visited_vns = (supermode == "LBP") ? zeros(Bool,N) : nothing

    if supermode == "LBP" || supermode == "RBP"       
        Ldn = zeros(N)
    else
        Ldn = nothing
    end

    phi = (mode == "TABL") ? lookupTable() : nothing

    Residues = (mode == "RBP" || mode == "RRBP") ? H*0.0 : nothing

    samples = (mode == "RRBP" && SAMPLESIZE != 0) ?
                Vector{Int}(undef,SAMPLESIZE) : nothing

    maxcoords, Factors, num_edges = (supermode == "RBP") ? 
        ([0,0], 1.0*H, sum(H)) : 
        (nothing,nothing,nothing)

################################## MAIN LOOP ###################################
    
    rng_noise = Xoshiro(rgn_seed_noise)

    rng_sample = (mode == "RRBP") ? Xoshiro(rng_seed_sample) : nothing

    ######################### FIRST RECEIVED SIGNAL ############################
    # In order to allow to test with a given received signal t_test, the first
    # received signal t must be set outside the main loop.
    if test
        if t_test === nothing # if no test signal was provided:
            # generate a new received signal
            received_signal!(t,noise,stdev,u,rng_noise)
        elseif length(t_test) != N
            # if a received test signal was given, but with wrong length
            throw(
                DimensionMismatch(
                    "length(t_test) should be $N, not $(length(t_test))"
                )
            )
        else
            # the received signal is the given test signal
            t = t_test
        end
    else
        # generate the first received signal outside the main loop
        received_signal!(t,noise,stdev,u,rng_noise)
    end        

    for j in 1:nreals

        # init the llr priors
        calc_Lf!(Lf,t,variance)
        if mode == "TABL"
            # scale for table
            Lf .*= SIZE_per_RANGE
        end            
        # initialize matrix Lr
        Lr .*= 0
        # initialize matrix Lq
        init_Lq!(Lq,Lf,vn2cn)

        if supermode == "RBP"
            # minsum_RBP_init!(Residues,maxcoords,Lq,signs,cn2vn)
            init_residues!(Residues,maxcoords,Lq,signs,cn2vn,Ms)
        end      
        # SPA routine
        # DECODED, i = BP!(
        DECODED = BP!(
                    supermode,
                    mode,
                    stop,
                    test,
                    max,
                    syndrome,
                    d,
                    c,
                    bit_error,
                    ber,
                    Lf,
                    Lq,
                    Lr,
                    cn2vn,
                    vn2cn,
                    Lrn,
                    signs,
                    phi,
                    printing,
                    Residues,
                    maxcoords,
                    Factors,
                    H,
                    rbpfactor,
                    num_edges,
                    Ldn,
                    visited_vns,
                    samples,
                    rng_sample,
                    Ms
                )                

        # bit error rate
        @fastmath @. BER += ber
        # iteration in which SPA stopped (iszero(syndrome) = true)
        # @inbounds iters[j] = i
        if !(DECODED)
            # frame error rate
            FER += 1
        end

        # received signal for the next realization (j+1)
        received_signal!(t,noise,stdev,u,rng_noise)

    end

    # @inbounds @fastmath FER /= nreals


    if test
        return Lr, Lq
    else
        # return log10.(FER), log10.(BER), iters
        return FER, BER
    end

end
