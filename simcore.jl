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
        msg::Vector{Bool},
        cword::Vector{Bool},
        snr::Real,
        H::BitMatrix,
        Zc::Integer,
        mode::String,
        supermode::String,
        flootype::String,
        alpha::AbstractFloat,
        trials::Integer,
        maxiter::Integer,
        stop::Bool,
        fast::Bool,
        rbpfactor::Union{AbstractFloat,Nothing},
        listsize1::Integer,
        listsize2::Integer,
        samplesize::Integer,
        rgn_seed_noise::Integer,
        rng_seed_sample::Integer,
        test::Bool,
        printtest::Bool
    )
    
################################## CONSTANTS ###################################
    
    # transform snr in standard deviations
    variance = 1 ./ (exp10.(snr/10))
    stdev = sqrt.(variance)
    
    u = Float64.(2*cword .- 1)
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
    syndrome = Vector{Bool}(undef,M)

    # prior llr (if mode == "MKAY" just the prior probabilities)
    Lf = (mode != "MKAY") ? zeros(N) : zeros(N,2)

    # noise
    L = length(cword)
    noise = Vector{Float64}(undef,L)

    # received signal
    signal = zeros(L)

    # bit-error
    biterror = Vector{Bool}(undef,N) 

    # Lq -> matrix of vn2cn messages (N x M)
    # Lr -> matrix of cn2vn messages (M x N)
    # if mode == "MKAY" the are matrices for bit = 0 and bit = 1
    Lq = (mode != "MKAY") ? zeros(N,M) : zeros(N,M,2)

    Lr = (mode != "MKAY") ? zeros(M,N) : zeros(M,N,2)

    Ms = (supermode == "RBP") ? H*0.0 : nothing
    
    # Set variables Lrn and signs depending on the floooding type (also used for dispatch)
    if (flootype == "TANH" && fast) || supermode == "LBP"
        Lrn = zeros(N)
        signs = nothing
    elseif flootype == "ALTN" || flootype == "TABL" 
        Lrn = zeros(N)
        signs = zeros(Bool,N)
    elseif flootype == "MSUM"
        Lrn = nothing
        signs = zeros(Bool,N)
    else
        Lrn = nothing
        signs = nothing
    end

   # Set other variables that depend on the mode

    visited_vns = (mode == "iLBP") ? zeros(Bool,N) : nothing

    if mode == "iLBP" || supermode == "RBP"       
        Ldn = zeros(N)
    else
        Ldn = nothing
    end

    phi = (mode == "TABL") ? lookupTable() : nothing

    Residues = (mode == "RBP" || mode == "Random-RBP") ? H*0.0 : nothing

    samples = (mode == "Random-RBP" && samplesize != 0) ?
                Vector{Int}(undef,samplesize) : nothing

    maxresidue = 0.0

    maxcoords, Factors, num_edges = (supermode == "RBP") ? 
        ([0,0], 1.0*H, sum(H)) : 
        (nothing,nothing,nothing)
    
    if mode == "List-RBP"
        listres1 = zeros(listsize1+1)
        listadd1 = zeros(Int,2,listsize1+1)
        listaddinv1 = zeros(Int,M,N)
        inlist1 = Matrix(false*H)
        if listsize2 != 0
            listres2 = zeros(listsize2+1)
            listadd2 = zeros(Int,2,listsize2+1)
        else
            listres2 = nothing
            listadd2 = nothing
        end
    else
        listres1 = nothing
        listadd1 = nothing
        listres2 = nothing
        listadd2 = nothing
        listaddinv1 = nothing
        inlist1 = nothing
        listsize2 = 0
    end 

################################## MAIN LOOP ###################################
    
    rng_noise = Xoshiro(rgn_seed_noise)

    rng_sample = (supermode == "RBP") ? Xoshiro(rng_seed_sample) : nothing  
    
    if Zc > 0
        cword = [msg[1:N-L]; cword]
    end

    for j in 1:trials

        if mode == "List-RBP"
            listres1 .*= 0.0
            listadd1 .*= 0
            @inbounds for m in 1:M
                for n in cn2vn[m]
                    inlist1[m,n] = false
                end
            end
        end

        # received signal for the next realization (j+1)
        received_signal!(signal,noise,stdev,u,rng_noise)

        bitvector .= true

        syndrome .= true

        decoded .= false
        
        # reinitialize matrices Lr and Lq
        @inbounds for m in 1:M
            for n in cn2vn[m]
                Lr[m,n] *= 0.0
                Lq[n,m] *= 0.0
            end
        end

        # init the llr priors
        calc_Lf!(view(Lf,2*Zc+1:N),signal,variance)
        if mode == "TABL"
            # scale for table
            Lf .*= SIZE_per_RANGE
        end
        init_Lq!(Lq,Lf,vn2cn)

        if supermode == "RBP"
            maxresidue = init_residues!(alpha,Residues,maxcoords,Lq,Lrn,signs,phi,cn2vn,Ms,
                listres1,listadd1,listaddinv1,listsize1,inlist1)
        end
            
        # SPA routine
        BP!(supermode,
            alpha,
            stop,
            test,
            maxiter,
            syndrome,
            bitvector,
            cword,
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
            listsize1,
            listsize2,
            listres1,
            listadd1,
            listres2,
            listadd2,
            listaddinv1,
            inlist1)                

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

