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
        codeword::Vector{Bool},
        snr::Real,
        H::BitMatrix,
        mode::String,
        supermode::String,
        trials::Integer,
        maxiter::Integer,
        stop::Bool,
        rbpfactor::Union{AbstractFloat,Nothing},
        rgn_seed_noise::Integer,
        rng_seed_sample::Integer,
        test::Bool,
        testsignal::Union{Vector{<:AbstractFloat},Nothing},
        printtest::Bool
    )
    
################################## CONSTANTS ###################################
    
    # transform snr in standard deviations
    variance = 1 ./ (exp10.(snr/10))
    stdev = sqrt.(variance)
    # BPKS
    u = Float64.(2*codeword .- 1)
    # list of checks and variables nodes
    vn2cn  = make_vn2cn_list(H)
    cn2vn  = make_cn2vn_list(H)

    ############################# PREALLOCATIONS ################################

    # frame error rate
    DECODED = zeros(Int,maxiter)
    decoded = Vector{Bool}(undef,maxiter)

    # bit error rate
    BER = zeros(Int,maxiter)
    ber = zeros(Int,maxiter)

    # estimate
    bitvector = zeros(Bool,N)

    # syndrome
    syndrome = ones(Bool,M)

    # prior llr (if mode == "MKAY" just the prior probabilities)
    Lf = (mode != "MKAY") ? Vector{Float64}(undef,N) : Matrix{Float64}(undef,N,2)

    # noise
    noise = Vector{Float64}(undef,N)

    # received signal
    signal = Vector{Float64}(undef,N)

    # bit-error
    biterror = Vector{Bool}(undef,N) 

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

    Residues = (mode == "RBP" || mode == "Random-RBP") ? H*0.0 : nothing

    samples = (mode == "Random-RBP" && SAMPLESIZE != 0) ?
                Vector{Int}(undef,SAMPLESIZE) : nothing

    maxcoords, Factors, num_edges = (supermode == "RBP") ? 
        ([0,0], 1.0*H, sum(H)) : 
        (nothing,nothing,nothing)

################################## MAIN LOOP ###################################
    
    rng_noise = Xoshiro(rgn_seed_noise)

    rng_sample = (mode == "Random-RBP") ? Xoshiro(rng_seed_sample) : nothing

    ######################### FIRST RECEIVED SIGNAL ############################
    # In order to allow to test with a given received signal, its first 
    # realization must be set outside the main loop.
    if test
        if testsignal === nothing # if no test signal was provided:
            # generate a new received signal
            received_signal!(signal,noise,stdev,u,rng_noise)
        elseif length(testsignal) != N
            # if a received test signal was given, but with wrong length
            throw(
                DimensionMismatch(
                    "length(testsignal) should be $N, not $(length(testsignal))"
                )
            )
        else
            # the received signal is the given test signal
            signal = testsignal
        end
    else
        # generate the first received signal outside the main loop
        received_signal!(signal,noise,stdev,u,rng_noise)
    end        

    for j in 1:trials

        # init the llr priors
        calc_Lf!(Lf,signal,variance)
        if mode == "TABL"
            # scale for table
            Lf .*= SIZE_per_RANGE
        end            
        # initialize matrix Lr
        Lr .*= 0
        # initialize matrix Lq
        init_Lq!(Lq,Lf,vn2cn)

        if mode == "RBP" || mode == "Random-RBP"
            init_residues!(Residues,maxcoords,Lq,signs,cn2vn,Ms)
            listres = nothing
            listadd = nothing
            inlist = nothing
        elseif mode == "Local-RBP"
            maxcoords = [1,cn2vn[1][1]]
            listres = nothing
            listadd = nothing
            inlist = nothing
        elseif mode == "List-RBP"
            listres = zeros(LISTSIZE)
            listres[1] = -1.0
            listadd = zeros(Int,2,LISTSIZE)
            listadd[1,1] = 1
            listadd[2,1] = cn2vn[1][1]
            inlist = Matrix(false*H)
            inlist[1,cn2vn[1][1]] = true       
        end
            
        # SPA routine
        decoded .= false
        BP!(supermode,
            mode,
            stop,
            test,
            maxiter,
            syndrome,
            bitvector,
            codeword,
            biterror,
            ber,
            decoded,
            Lf,
            Lq,
            Lr,
            Ms,
            cn2vn,
            vn2cn,
            Lrn,
            signs,
            phi,
            printtest,
            Residues,
            maxcoords,
            Factors,
            rbpfactor,
            num_edges,
            Ldn,
            visited_vns,
            samples,
            rng_sample,
            listres,
            listadd,
            LISTSIZE,
            inlist)                

        # bit error rate
        @. BER += ber
        @. DECODED += decoded

        # received signal for the next realization (j+1)
        received_signal!(signal,noise,stdev,u,rng_noise)

    end

    if test
        return Lr, Lq
    else
        return DECODED, BER
    end

end

