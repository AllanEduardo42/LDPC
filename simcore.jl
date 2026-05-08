################################################################################
# Allan Eduardo Feitosa
# 25 Feb 2025
# Core routine to estimate the LPCD performance (FER sum_ber x SNR)

# simcore functions
include("./Simulation core functions/auxiliary_functions.jl")
include("./Simulation core functions/encode_LDPC.jl")
include("./Simulation core functions/tanh_V2C.jl")
include("./Simulation core functions/calc_syndrome.jl")
include("./Simulation core functions/calc_C2V.jl")

function simcore(
    A::Int,                         # payload length
    K::Int,                         # payload + CRC length
    R::Float64,                     # effective rate
    G::Int,                         # transmitted signal length
    g_CRC::Vector{Bool},            # CRC polynomial
    ebn0::Float64,                  # EbN0
    H::Matrix{Bool},                # Parity-Check Matrix
    H1::Matrix{Bool},               # H = [H1 H2], size(H1) = (M,K)
    L::Matrix{Bool},                # L*U = H2  (PEG)
    U::Matrix{Bool},                # L*U = H2  (PEG)
    Nc::Vector{Vector{Int}},        # Nc[ci] : neighborhood of check node ci
    Nv::Vector{Vector{Int}},        # Nv[vj] : neighborhood of variable node vj
    eH::Matrix{Int},                # Exponential Matrix (5GNR and WiMAX)
    protocol::String,               # PEG, 5GNR or WiMAX
    LS::Int,                        # 5GNR and WiMAX matrix liftsize
    algorithm::String,              # BP algorithm (flooding, RBP etc.)
    bptype::String,                 # Type of BP implementation (FAST, TANH etc.)
    max_frame_errors::Int,                # Number of max frame errors
    maxiter::Int,                   # Maximum number of BP iterations
    rayleigh::Bool,                 # Rayleigh fading flag
    C_DR_iter::Int,                 # C&DR switch iter
    decayfactor::Float64,           # RBP decay factor
    listsizes::Vector{Int},         # sizes of the list for List-RBP algorithm
    ci_gamma::Float64,              # CI gamma threshold
    rgn_seed::Int,                  # random seed to generate noise and message
    test::Bool,                     # if "true", perform test algorithm
    printtest::Bool                 # if "true", print test algorithm results
)
    
################################## CONSTANTS ###################################

    # Parity-Check Matrix dimensions
    M,N = size(H)

    # 2*LS in NR 5G (determines the number of initial punctured bits)
    # = 0 otherwise
    twoLs = 0

    # Parity Check Matrix constants
    if protocol == "WiMAX" || protocol == "5GNR"
        eM, eN = size(eH)
        eK = eN - eM 
        if protocol == "WiMAX"            
            eB = eM  
            S = 0   
        elseif protocol == "5GNR"
            twoLs = 2*LS
            eB = 4
            S = eK*LS - K               # filler bits
        end
    end

    # number of edges in the graph
    num_edges = sum(H)

    # transform EbN0 in standard deviations
    variance = exp10.(-ebn0/10) / (2*R)
    stdev = sqrt.(variance)

    # Set the random seeds
    rgn = Xoshiro(rgn_seed)

    # min-sum constants
    msum2 = false
    if bptype == "MSUM"
        msum_factor = ALPHA
    elseif bptype == "MSUMRBP"
        msum2 = true
        msum_factor = ALPHA2
    else
        msum_factor = nothing
    end

################################# PREALLOCATIONS ###############################

    msg = Vector{Bool}(undef,A)         # payload
    
    b = Vector{Bool}(undef,K)           # payload + CRC

    cword = Vector{Bool}(undef,N)       # codeword (payload + CRC + parity bits)

    if protocol == "PEG"
        Cw = cword
        w = Vector{Bool}(undef,M)       # parity bits
    else
        Cw = Matrix{Bool}(undef,LS,eN)  # payload + CRC + filler bits + parity bits
        aux1 = Vector{Bool}(undef,LS)
        aux2 = Vector{Bool}(undef,LS)
        aux3 = Vector{Bool}(undef,LS)    
        W = Matrix{Bool}(undef,LS,eB)    
    end

    # frame error rate
    sum_decoded = zeros(Int,maxiter)
    decoded = Vector{Bool}(undef,maxiter)

    # bit error rate
    sum_ber = zeros(Int,maxiter)
    ber = Vector{Int}(undef,maxiter)

    # estimate
    bitvector = Vector{Bool}(undef,N)    

    # prior llr
    prior_LLRs = Vector{Float64}(undef,N)

    # if there are punctured bits
    if twoLs > 0
        for i in 1:twoLs
            prior_LLRs[i] = 0.0
        end
    end

    # V2C -> matrix of V2C messages (N x M)
    # C2V -> matrix of C2V messages (N x M)
    V2C = Matrix{Float64}(undef,M,N)
    C2V = Matrix{Float64}(undef,M,N)

    # received signal
    signal = Vector{Float64}(undef,G)

    # bit-error
    biterror = Vector{Bool}(undef,N) 

    # Rayleigh fading channel
    if rayleigh
        x1 = Vector{Float64}(undef,G)
        x2 = Vector{Float64}(undef,G)
        fading = Vector{Float64}(undef,G)   # fading .= sqrt(x_1^2 + x_2^2)
    else
        x1 = nothing
        x2 = nothing
        fading = nothing
    end

    if test
        syndrome = Vector{Bool}(undef,M)
    end

############################### RBP PREALLOCATIONS #############################
    
    # RBP switchs 
    RBP_based = algorithm != "Flooding" && algorithm != "LBP"      

    if RBP_based

        # Consensus
        RBP_consensus = algorithm == "C-RBP"    || 
                        algorithm == "C&R-RBP"  ||
                        algorithm == "C&DR-RBP" 

        # Algorithms that used the RBP residual decay strategy
        RBP_decay = RBP_consensus           ||
                    algorithm == "RD-RBP"   || 
                    algorithm == "List-RBP"
                     

        # Algorithms that uses the default RBP initialization
        RBP_init = algorithm != "List-RBP" &&
                   algorithm != "OV-RBP"   &&
                   algorithm != "VC-RBP"
        
        # newly calculate C2V messages
        newC2V = Matrix{Float64}(undef,M,N)

        if algorithm != "OV-RBP"
            Residuals = Matrix{Float64}(undef,M,N)
        else
            Residuals = Vector{Float64}(undef,N)
        end

        if RBP_decay
            Factors = Matrix{Float64}(undef,M,N)
        else
            Factors = nothing
        end

        if RBP_init
            alpha = Vector{Float64}(undef,M)
        end

        # default number of loop repetitions
        num_reps = num_edges

        if RBP_consensus
            switch_R = false            # switch Return
            switch_C_DR = false         # switch Consensus & Return
            if algorithm == "C&R-RBP"
                num_reps -= N
                switch_R = true
            elseif algorithm == "C&DR-RBP"
                switch_C_DR = true
            end
        elseif algorithm == "NW-RBP"
            num_reps = M  
        elseif algorithm == "List-RBP"
            listsize = listsizes[1]
            list = Vector{Float64}(undef,listsize+1)
            list[end] = 0.0
            coords = Matrix{Int}(undef,3,listsize+1)
            inlist = zeros(Bool,M,N)
            listsize2 = listsizes[2]
            local_list = Vector{Float64}(undef,listsize2+1)
            local_coords = Matrix{Int}(undef,3,listsize2+1)
        elseif algorithm == "VC-RBP"
            alpha = Vector{Float64}(undef,N) 
        elseif algorithm == "OV-RBP"
            num_reps = N   
            LLRs = Vector{Float64}(undef,N)
            newLLRs = Vector{Float64}(undef,N)
            C = zeros(Bool,N)
            C1 = zeros(Bool,N)
            upc = zeros(Int,N)
        elseif algorithm == "CI-RBP"
            Prob = Vector{Float64}(undef,N)     # Probability of bj = 0
            Dn = Vector{Float64}(undef,N)
        elseif algorithm == "UBP-RBP"
            UBP = ones(Bool,M)
        elseif algorithm == "RBP-D1VN"
            F = zeros(Bool,N)
            degree_vn = sum(H,dims=1)[:]
            E1 = 0
            for dn in degree_vn
                if dn == 1
                    E1 += 1
                end
            end
            num_reps -= E1
        end
    end

################################## MAIN LOOP ###################################
    
    frame_errors = 0
    trials = 0

    #while frame_errors < max_frame_errors
    @fastmath @inbounds while frame_errors < max_frame_errors

        trials += 1

        ### 1) generate the random message
        rand!(rgn,msg,Bool)

        ### 2) generate the codeword
        append_CRC!(Cw,b,msg,g_CRC,A,K)
        if protocol == "PEG"
            encode_LDPC_LU!(Cw,H1,w,L,U,M,K)
            for i in 1:M
                Cw[K+i] = w[i]
            end
        else
            for i in K+1:K+S
                Cw[i] = false           # filler bits (they are removed later)
            end
            encode_LDPC_BG!(cword,Cw,W,aux1,aux2,aux3,K,N,eH,eM,eK,eB,S,LS)
        end
        # if in test algorithm, verify the encoding
        if test
            _gf2_mat_mult!(syndrome,H,cword,M,N)
            if !iszero(syndrome)
                throw(error("encoding error"))
            end
        end

        ### 3) sum the noise to the modulated cword to produce the received signal
        received_signal!(signal,cword,G,twoLs,stdev,rgn,rayleigh,fading,x1,x2)
        
        ### 4) reset simulation variables
        decoded .= false
        resetmatrix!(C2V,Nv,0.0)

        ### 5) init the LLR priors
        calc_prior_LLRs!(prior_LLRs,twoLs,signal,variance,rayleigh,fading,bptype)

        # initial estimate of the bits (0-th iteration)
        for i in eachindex(bitvector)
            bitvector[i] = signbit(prior_LLRs[i])
        end

        ### 6) init the V2C matrix
        init_V2C!(V2C,prior_LLRs,Nv,msum_factor)

        ### 6.5) print test info
        if printtest
            println("""
________________________________________________________________________________

                                    TRIAL #$trials")
________________________________________________________________________________
""")
            print_test("Message",msg)
            print_test("Transmitted codeword",cword[twoLs+1:end])
            println()
            println("                              ### Iteration #0 ###")
            calc_syndrome!(syndrome,bitvector,Nc)
            biterror .= (bitvector .≠ cword)
            print_test("Syndrome",syndrome)  
            println("Syndrome rate: $(sum(syndrome))/$M")
            print_test("Bit error",biterror)   
            println("Bit error rate: $(sum(biterror))/$N")            
            println()
        end

        ### 7) init the RBP methods
        if RBP_based

            if RBP_decay
                resetmatrix!(Factors,Nv,1.0)
            end

            if RBP_init
                init_residuals!(V2C,Nc,newC2V,alpha,Residuals,msum_factor)
                if algorithm == "CI-RBP"
                    for vj in eachindex(Nv)
                        Prob[vj] = calc_prob(prior_LLRs[vj])
                        calc_Dn!(Dn,Prob,newC2V,prior_LLRs,vj,Nv)
                    end
                end
            elseif algorithm == "List-RBP"
                list .= 0.0
                coords .= 0
                resetmatrix!(inlist,Nv,false)
                local_list .= 0.0
                local_coords .= 0
                init_list!(V2C,Nc,msum_factor,newC2V,inlist,Residuals,list,
                                                                coords,listsize)  
            elseif algorithm == "VC-RBP"
                init_VC!(V2C,Nv,alpha,Residuals)
            elseif algorithm == "OV-RBP"
                max_upc = init_OV!(V2C,Nc,Nv,newC2V,Residuals,msum_factor,
                prior_LLRs,LLRs,newLLRs,C,C1,upc)
            end
        end

        ### 8) BP routine

        rbp_not_converged = true        # avoids useless iterations if false

        iter = 0
        while iter < maxiter

            iter += 1
    
            if RBP_based
                if algorithm == "RBP"
                    rbp_not_converged = RBP!(
                        bitvector,
                        V2C,
                        C2V,
                        prior_LLRs,
                        Nc,
                        Nv,
                        msum_factor,
                        msum2,
                        num_reps,
                        newC2V,
                        Residuals,
                        alpha
                        )
                elseif algorithm == "RD-RBP"
                    rbp_not_converged = RD_RBP!(
                        bitvector,
                        V2C,
                        C2V,
                        prior_LLRs,
                        Nc,
                        Nv,
                        msum_factor,
                        msum2,
                        num_reps,
                        newC2V,
                        Residuals,
                        alpha,
                        decayfactor,
                        Factors
                        )
                elseif RBP_consensus
                    if switch_C_DR && iter ≥ C_DR_iter
                        switch_R = true
                        switch_C_DR = false
                        num_reps -= N
                    end             
                    rbp_not_converged = C_and_R_RBP!(
                        bitvector,
                        V2C,
                        C2V,
                        prior_LLRs,
                        Nc,
                        Nv,
                        msum_factor,
                        msum2,
                        num_reps,
                        newC2V,
                        Residuals,
                        alpha,
                        decayfactor,
                        Factors,
                        switch_R
                        )
                elseif algorithm == "NW-RBP"
                    rbp_not_converged = NW_RBP!(
                        bitvector,
                        V2C,
                        C2V,
                        prior_LLRs,
                        Nc,
                        Nv,
                        msum_factor,
                        msum2,
                        num_reps,
                        newC2V,
                        alpha                    
                    )
                elseif algorithm == "SVNF"
                    rbp_not_converged = SVNF!(
                        bitvector,
                        V2C,
                        C2V,
                        prior_LLRs,
                        Nc,
                        Nv,
                        msum_factor,
                        msum2,
                        num_reps,
                        newC2V,
                        Residuals,
                        twoLs,
                        N
                    )
                elseif algorithm == "List-RBP"
                    rbp_not_converged = List_RBP!(
                        bitvector,
                        V2C,
                        C2V,
                        prior_LLRs,
                        Nc,
                        Nv,
                        msum_factor,
                        msum2,
                        num_reps,
                        newC2V,
                        Residuals,
                        decayfactor,
                        Factors,
                        coords,
                        inlist,
                        list,
                        local_list,
                        local_coords,
                        listsize,
                        listsize2
                    )
                elseif algorithm == "VC-RBP"
                    rbp_not_converged = VC_RBP!(
                        bitvector,
                        V2C,
                        C2V,
                        prior_LLRs,
                        Nc,
                        Nv,
                        msum_factor,
                        msum2,
                        num_reps,
                        Residuals,
                        alpha
                    )
                elseif algorithm == "OV-RBP"
                    rbp_not_converged = OV_RBP!(
                        bitvector,
                        V2C,
                        C2V,
                        prior_LLRs,
                        Nc,
                        Nv,
                        msum_factor,
                        msum2,
                        num_reps,
                        newC2V,
                        Residuals,                    
                        LLRs,
                        newLLRs,
                        C,
                        C1,
                        upc,
                        max_upc
                    )
                elseif algorithm == "CI-RBP"
                    rbp_not_converged = CI_RBP!(
                        bitvector,
                        V2C,
                        C2V,
                        prior_LLRs,
                        Nc,
                        Nv,
                        msum_factor,
                        msum2,
                        num_reps,
                        newC2V,
                        Residuals,
                        alpha,
                        Dn,
                        Prob,
                        ci_gamma
                        )
                elseif algorithm == "UBP-RBP"
                    rbp_not_converged = UBP_RBP!(
                        bitvector,
                        V2C,
                        C2V,
                        prior_LLRs,
                        Nc,
                        Nv,
                        msum_factor,
                        msum2,
                        num_reps,
                        newC2V,
                        Residuals,
                        alpha,
                        UBP
                        )
                    elseif algorithm == "RBP-D1VN"
                        rbp_not_converged = RBP_D1VN!(
                        bitvector,
                        V2C,
                        C2V,
                        prior_LLRs,
                        Nc,
                        Nv,
                        msum_factor,
                        msum2,
                        num_reps,
                        newC2V,
                        Residuals,
                        alpha,
                        F,
                        degree_vn
                        )
                else
                    throw(error(lazy"Invalid RBP-based Algorithm"))  
                end

                # Reset Factors
                if RBP_decay
                    resetmatrix!(Factors,Nv,1.0)
                end

            elseif algorithm == "Flooding"
                flooding!(
                    bitvector,
                    V2C,
                    C2V,
                    prior_LLRs,
                    Nc,
                    Nv,
                    msum_factor
                    )
            elseif algorithm == "LBP"
                LBP!(
                    bitvector,
                    V2C,
                    C2V,
                    prior_LLRs,
                    Nc,
                    Nv,
                    msum_factor
                    )
            end            

            # Evalute the bit error vector
            for vj in 1:N
                biterror[vj] = bitvector[vj] ≠ cword[vj]
            end
            # Evaluate BER
            ber[iter] = sum(biterror)

            # print test info
            if printtest
                println("                              ### Iteration #$iter ###")       
                calc_syndrome!(syndrome,bitvector,Nc) 
                print_test("Syndrome",syndrome)  
                println("Syndrome rate: $(sum(syndrome))/$M")
                print_test("Bit error",biterror)   
                println("Bit error rate: $(sum(biterror))/$N")
                println()
            end

            # Verify if all check equations are satisfied (syndrome == zero)
            if iszerosyndrome(bitvector,Nc)
                if iszero(biterror)
                    decoded[iter] = true
                end
                if printtest
                    println("#### Algorithm stopped: zero syndrome at iteration $iter ####")
                    println()
                end
                break                   # stops the algorithm
            end
            if !rbp_not_converged       # all residuals are zero
                if printtest
                    println("#### Algorithm stopped: BP converged at iteration $iter ####")
                    println()
                end
                break
            end
        end

        if iter < maxiter               # if the algorithm stopped before maxiter
            for i = iter+1:maxiter
                decoded[i] = decoded[iter]
                ber[i] = ber[iter]
            end
        end

        # reset C&DR-RBP
        if algorithm == "C&DR-RBP"
            switch_R = false
            switch_C_DR = true
            num_reps = num_edges
        end

        ### 9) bit error rate
        for i=1:maxiter
            sum_ber[i] += ber[i]
            sum_decoded[i] += decoded[i]
        end

        if !decoded[maxiter]
            frame_errors += 1
        end

        if printtest
            println("Frame Errors: $(frame_errors)/$(max_frame_errors)")
            println()
        end
    end

    if printtest
        println("Monte Carlo trials ended: maximum number of frame errors reached.")
        println("Total number of trials: $trials")
    end

    if test
        return C2V, V2C, sum_decoded, sum_ber, trials
    else
        return nothing, nothing, sum_decoded, sum_ber, trials
    end

end