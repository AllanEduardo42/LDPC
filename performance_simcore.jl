################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Core routine to estimate the LPCD performance (FER BER x SNR)

include("auxiliary_functions.jl")
include("lookupTable.jl")
include("BP.jl")
include("calc_Lf.jl")

function
    performance_simcore(
        message::Vector{Bool},
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

    ############################# PREALLOCATIONS ###############################

    # frame error rate
    DECODED = zeros(Int,maxiter)
    decoded = Vector{Bool}(undef,maxiter)

    # bit error rate
    BER = zeros(Int,maxiter)
    ber = zeros(Int,maxiter)

    # estimate
    bitvector = Vector{Bool}(undef,N)

    # syndrome
    syndrome = Vector{Bool}(undef,N)

    # prior llr (if mode == "MKAY" just the prior probabilities)
    Lf = (mode != "MKAY") ? zeros(N) : zeros(N,2)

    # noise
    L = length(codeword)
    noise = Vector{Float64}(undef,L)

    # received signal
    signal = zeros(L)
    NmL = N - L
    if NmL > 0
        Lf[1:NmL] .= -2*eps()/variance
    end

    # bit-error
    biterror = Vector{Bool}(undef,N) 

    # Lq -> matrix of vn2cn messages (N x M)
    # Lr -> matrix of cn2vn messages (M x N)
    # if mode == "MKAY" the are matrices for bit = 0 and bit = 1
    Lq = (mode != "MKAY") ? zeros(N,M) : zeros(N,M,2)

    Lr = (mode != "MKAY") ? zeros(M,N) : zeros(M,N,2)

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

    maxresidue = 0.0

    maxcoords, Factors, num_edges = (supermode == "RBP") ? 
        ([0,0], 1.0*H, sum(H)) : 
        (nothing,nothing,nothing)
    
    if mode == "List-RBP"
        listres = zeros(LISTSIZE+1)
        listadd = zeros(Int,2,LISTSIZE+1)
        inlist = Matrix(false*H)
    else
        listres = nothing
        listadd = nothing
        inlist = nothing
    end 

################################## MAIN LOOP ###################################
    
    rng_noise = Xoshiro(rgn_seed_noise)

    rng_sample = (mode == "Random-RBP") ? Xoshiro(rng_seed_sample) : nothing  
    
    if NmL > 0
        codeword = [message[1:N-L]; codeword]
    end

    for j in 1:trials

        # received signal for the next realization (j+1)
        received_signal!(signal,noise,stdev,u,rng_noise)

        bitvector .= true

        syndrome .= true

        decoded .= false
        
        # reinitialize matrix Lr
        Lr *= 0.0

        # init the llr priors
        calc_Lf!(view(Lf,NmL+1:N),signal,variance)
        if mode == "TABL"
            # scale for table
            Lf .*= SIZE_per_RANGE
        end
        # reinitialize matrix Lq
        Lq *= 0.0
        init_Lq!(Lq,Lf,vn2cn)

        if supermode == "RBP"
            maxresidue = init_residues!(Residues,maxcoords,Lq,signs,cn2vn,Ms,
                listres,listadd,LISTSIZE,inlist)   
        end
            
        # SPA routine
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
            maxresidue,
            maxcoords,
            Factors,
            rbpfactor,
            num_edges,
            Ldn,
            visited_vns,
            samples,
            rng_sample,
            LISTSIZE,
            listres,
            listadd,
            inlist)                

        # bit error rate
        @. BER += ber
        @. DECODED += decoded

    end

    if test
        return Lr, Lq
    else
        return DECODED, BER
    end

end

