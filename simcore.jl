################################################################################
# Allan Eduardo Feitosa
# 25 Feb 2025
# Core routine to estimate the LPCD performance (FER sum_ber x SNR)

include("auxiliary_functions.jl")
include("lookupTable.jl")
include("calc_prior_LLRs.jl")
include("calc_syndrome.jl")
include("flooding.jl")
include("LBP.jl")
include("RBP.jl")
include("RD-RBP.jl")
include("C&R-RBP.jl")
include("VC-RBP.jl")
include("OV-RBP.jl")
include("E_NV_RBP.jl")
include("List-RBP.jl")
include("SVNF.jl")
include("NW-RBP.jl")
include("CI-RBP.jl")
include("UBP-RBP.jl")
include("RBP_D1VN.jl")
include("./List functions/init_list.jl")
include("./RBP functions/init_RBP.jl")
include("update_C2V.jl")
include("update_V2C.jl")
include("NR LDPC/NR_LDPC_functions.jl")
include("encode_LDPC.jl")
include("tanh_V2C.jl")

function
    simcore(
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
        eH::Matrix{Int},                # Exponential Matrix (NR5G and WiMAX)
        protocol::String,               # PEG, NR5G or WiMAX
        LS::Int,                        # NR5G and WiMAX matrix liftsize
        algorithm::String,              # BP algorithm (flooding, RBP etc.)
        bptype::String,                 # Type of BP implementation (FAST, TANH etc.)
        trials::Int,                    # Number of trials
        maxiter::Int,                   # Maximum number of BP iterations
        stop::Bool,                     # "true" if routine stops at zero syndrome
        rayleigh::Bool,
        C_DR_iter::Int,                 # C&DR switch iter
        decayfactor::Float64,           # RBP decay factor
        listsizes::Vector{Int},         # sizes of the list for List-RBP algorithm
        ci_gamma::Float64,              # CI gamma threshold
        rgn_seed::Int,                  # random seed to generate noise and message
        test::Bool,                     # if "true", perform test algorithm
        printtest::Bool                # if "true", print test algorithm results
    )::Tuple{Union{Matrix{Float64},Nothing},Union{Matrix{Float64},Nothing},Vector{Int},Vector{Int}}
    
################################## CONSTANTS ###################################

    # Parity-Check Matrix dimensions
    M,N = size(H)

    # 2*LS in NR 5G (determines the number of initial punctured bits)
    # = 0 otherwise
    twoLs = 0

    # Parity Check Matrix constants
    if protocol == "WiMAX" || protocol == "NR5G"
        eM, eN = size(eH)
        eK = eN - eM 
        if protocol == "WiMAX"            
            eB = eM  
            S = 0   
        elseif protocol == "NR5G"
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

    # prior llr (if algorithm == "MKAY" it is just the prior probabilities)
    prior_LLRs = (bptype != "MKAY") ? Vector{Float64}(undef,N) : Matrix{Float64}(undef,N,2)

    if twoLs > 0
        for i in 1:twoLs
            prior_LLRs[i] = 0.0
        end
    end

    # received signal
    signal = Vector{Float64}(undef,G)

    # bit-error
    biterror = Vector{Bool}(undef,N) 

    if rayleigh
        x1 = Vector{Float64}(undef,G)
        x2 = Vector{Float64}(undef,G)
        fading = Vector{Float64}(undef,G)   # fading .= sqrt(x_1^2 + x_2^2)
    else
        x1 = nothing
        x2 = nothing
        fading = nothing
    end

    # V2C -> matrix of V2C messages (N x M)
    # C2V -> matrix of C2V messages (N x M)
    # if algorithm == "MKAY" the are different matrices for bit = 0 and bit = 1
    V2C = (bptype != "MKAY") ? Matrix{Float64}(undef,M,N) : Array{Float64,3}(undef,M,N,2)
    C2V = (bptype != "MKAY") ? Matrix{Float64}(undef,M,N) : Array{Float64,3}(undef,M,N,2)
    
    msum2 = false
    if bptype == "MSUM"
        msum_factor = ALPHA
    elseif bptype == "MSUM2"
        msum2 = true
        msum_factor = ALPHA2
    else
        msum_factor = nothing
    end

    phi = (bptype == "TABL") ? lookupTable() : nothing

    if test
        syndrome = Vector{Bool}(undef,M)
    end

############################### RBP PREALLOCATIONS #############################
    
    # RBP switchs 
    RBP_based = algorithm != "Flooding" && 
                algorithm != "LBP"      

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
        RBP_init = RBP_consensus          || 
                   algorithm == "RBP"     ||
                   algorithm == "RD-RBP"  ||
                   algorithm == "NW-RBP"  || 
                   algorithm == "SVNF"    ||
                   algorithm == "CI-RBP"  ||
                   algorithm == "UBP-RBP" ||
                   algorithm == "RBP-D1VN"

        
        # newly calculate C2V messages
        newC2V = Matrix{Float64}(undef,M,N)

        if algorithm != "OV-RBP"
            Residues = Matrix{Float64}(undef,M,N)
        else
            Residues = Vector{Float64}(undef,N)
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
            Prob0 = Vector{Float64}(undef,N)
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
    
    # for trial in 1:trials
    @fastmath @inbounds for trial in 1:trials

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

        # print info if in test algorithm
        if test && printtest
            println("Trial #$trial:")
            print_test("Message",msg)
            print_test("Transmitted codeword",cword[twoLs+1:end])
        end
        
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
        # print info if in test algorithm
        if test && printtest
            println()
            println("### Iteration #0 ###")
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
                init_residues!(V2C,Nc,phi,newC2V,alpha,Residues,msum_factor)
                if algorithm == "CI-RBP"
                    for vj in eachindex(Nv)
                        Prob0[vj] = calc_prob(prior_LLRs[vj])
                        calc_Dn!(Dn,Prob0,newC2V,prior_LLRs,vj,Nv)
                    end
                end
            elseif algorithm == "List-RBP"
                list .= 0.0
                coords .= 0
                resetmatrix!(inlist,Nv,false)
                local_list .= 0.0
                local_coords .= 0
                init_list!(V2C,Nc,phi,msum_factor,newC2V,inlist,
                                                Residues,list,coords,listsize)  
            elseif algorithm == "VC-RBP"
                for vj in eachindex(Nv)
                    alp = 0.0
                    for ci in Nv[vj]
                        li = LinearIndices(V2C)[ci,vj]
                        residue = abs(V2C[li])
                        if residue > alp
                            alp = residue
                        end
                        Residues[li] = residue
                    end
                    alpha[vj] = alp
                end
            elseif algorithm == "OV-RBP"
                LLRs = copy(prior_LLRs)
                # 3 - 7
                for ci in eachindex(Nc)
                    Nci = Nc[ci]
                    for vj in Nci
                        newC2V[ci,vj] = calc_C2V(Nci,ci,vj,V2C,msum_factor) 
                    end
                end
                #8 - 10                    
                C .= false                      # set C of VNs whose sign of the LLR changes
                upc .= 0                        # unsatisfied parity check equations
                max_upc = 0                     # maximum number of unsatified parity check equations 
                for vj in eachindex(Nv)
                    oldllr = LLRs[vj]
                    newllr = calc_post_LLR(vj,Nv[vj],prior_LLRs,newC2V)
                    newLLRs[vj] = newllr
                    Residues[vj] = abs(newllr - oldllr)        
                    if sign(oldllr)*sign(newllr) < 0
                        C[vj] = true
                        count_upc = 0
                        for ci in Nv[vj]
                            if _calc_syndrome(bitvector,Nc[ci])
                                count_upc += 1
                            end
                        end
                        upc[vj] = count_upc
                        if count_upc > max_upc
                            max_upc = count_upc
                        end
                    end
                end

                # 11
                C1 .= false                     # set of VNs with upc = max_upc
                for vj in eachindex(Nv)
                    p = upc[vj]
                    if C[vj] && (p == max_upc)  # every VN in C1 is in C
                        C1[vj] = true
                    end
                end       
            end
        end

        ### 8) BP routine
        zero_syn = false
        rbp_not_converged = true

        iter = 0
        while iter < maxiter && rbp_not_converged && !zero_syn

            iter += 1

            if test && printtest  
                println("### Iteration #$iter ###")
            end
    
            if RBP_based
                if algorithm == "RBP"
                    rbp_not_converged = RBP!(
                        bitvector,
                        V2C,
                        C2V,
                        prior_LLRs,
                        Nc,
                        Nv,
                        phi,
                        msum_factor,
                        msum2,
                        num_reps,
                        newC2V,
                        Residues,
                        alpha,
                        rbp_not_converged
                        )
                elseif algorithm == "RD-RBP"
                    rbp_not_converged = RD_RBP!(
                        bitvector,
                        V2C,
                        C2V,
                        prior_LLRs,
                        Nc,
                        Nv,
                        phi,
                        msum_factor,
                        msum2,
                        num_reps,
                        newC2V,
                        Residues,
                        alpha,
                        decayfactor,
                        Factors,
                        rbp_not_converged
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
                        phi,
                        msum_factor,
                        msum2,
                        num_reps,
                        newC2V,
                        Residues,
                        alpha,
                        decayfactor,
                        Factors,
                        rbp_not_converged,
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
                        phi,
                        msum_factor,
                        msum2,
                        num_reps,
                        newC2V,
                        alpha,
                        rbp_not_converged                    
                    )
                elseif algorithm == "SVNF"
                    rbp_not_converged = SVNF!(
                        bitvector,
                        V2C,
                        C2V,
                        prior_LLRs,
                        Nc,
                        Nv,
                        phi,
                        msum_factor,
                        msum2,
                        num_reps,
                        newC2V,
                        Residues,
                        rbp_not_converged,
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
                        phi,
                        msum_factor,
                        msum2,
                        num_reps,
                        newC2V,
                        Residues,
                        decayfactor,
                        Factors,
                        rbp_not_converged,
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
                        phi,
                        msum_factor,
                        msum2,
                        num_reps,
                        Residues,
                        alpha,                    
                        rbp_not_converged
                    )
                elseif algorithm == "OV-RBP"
                    rbp_not_converged = OV_RBP!(
                        bitvector,
                        V2C,
                        C2V,
                        prior_LLRs,
                        Nc,
                        Nv,
                        phi,
                        msum_factor,
                        msum2,
                        num_reps,
                        newC2V,
                        Residues,
                        rbp_not_converged,                    
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
                        phi,
                        msum_factor,
                        msum2,
                        num_reps,
                        newC2V,
                        Residues,
                        alpha,
                        rbp_not_converged,
                        Dn,
                        Prob0,
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
                        phi,
                        msum_factor,
                        msum2,
                        num_reps,
                        newC2V,
                        Residues,
                        alpha,
                        rbp_not_converged,
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
                            phi,
                            msum_factor,
                            msum2,
                            num_reps,
                            newC2V,
                            Residues,
                            alpha,
                            rbp_not_converged,
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
                    phi,
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
                    phi,
                    msum_factor
                    )
            end            

            # Evalute the bit error vector
            for vj in 1:N
                biterror[vj] = bitvector[vj] ≠ cword[vj]
            end
            # Evaluate BER
            ber[iter] = sum(biterror)

            # print info if in test algorithm
            if test && printtest        
                calc_syndrome!(syndrome,bitvector,Nc) 
                print_test("Syndrome",syndrome)  
                println("Syndrome rate: $(sum(syndrome))/$M")
                print_test("Bit error",biterror)   
                println("Bit error rate: $(sum(biterror))/$N")
                println()
            end
            if stop
            # Verify if all check equations are satisfied (syndrome vector is zero)
                zero_syn = iszerosyndrome(bitvector,Nc)
                if zero_syn || !rbp_not_converged
                    if iszero(biterror)
                        decoded[iter] = true
                    end
                    if test && printtest
                        if zero_syn
                            println("#### Zero Syndrome at iteration $iter ####")
                        end
                        if !rbp_not_converged
                            println("#### BP converged at iteration $iter ####")
                        end
                    end
                end
            else
                if iszero(biterror)
                    decoded[iter] = true
                end
            end

        end

        if iter < maxiter
            for i = iter+1:maxiter
                decoded[i] = decoded[iter]
                ber[i] = ber[iter]
            end
        end

        # reset for C&DR-RBP
        if algorithm == "C&DR-RBP"
            switch_R = false
            switch_C_DR = true
            num_reps = num_edges
        end

        ###  bit error rate
        for i=1:maxiter
            sum_ber[i] += ber[i]
            sum_decoded[i] += decoded[i]
        end
    end

    # MKAY compatibility
    if test && bptype == "MKAY"
        retr = Matrix{Float64}(undef,M,N)
        retq = Matrix{Float64}(undef,M,N)
        for ci in eachindex(Nc)
            for vj in Nc[ci]
                retr[ci,vj] = log.(C2V[ci,vj,1]) - log.(C2V[ci,vj,2])
                retq[ci,vj] = log.(V2C[ci,vj,1]) - log.(V2C[ci,vj,2])
            end
        end
        C2V = retr
        V2C = retq
    end

    if test
        return C2V, V2C, sum_decoded, sum_ber
    else
        return nothing, nothing, sum_decoded, sum_ber
    end

end