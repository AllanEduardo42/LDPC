################################################################################
# Allan Eduardo Feitosa
# 25 Feb 2025
# Core routine to estimate the LPCD performance (FER sum_ber x SNR)

include("auxiliary_functions.jl")
include("lookupTable.jl")
include("calc_Lf.jl")
include("calc_syndrome.jl")
include("flooding.jl")
include("LBP.jl")
include("RBP.jl")
include("VC-RBP.jl")
include("OV-RBP.jl")
include("E_NV_RBP.jl")
include("List_RBP.jl")
include("SVNF.jl")
include("D-SVNF.jl")
include("NW_RBP.jl")
include("RPD.jl")
include("./List functions/init_list.jl")
include("./RBP functions/init_RBP.jl")
include("./RBP functions/init_SVNF.jl")
include("./RBP functions/init_NW_RBP.jl")
include("./RBP functions/init_LD_RBP.jl")
include("update_Lq.jl")
include("update_Lr.jl")
include("NR LDPC/NR_LDPC_functions.jl")
include("encode_LDPC.jl")
include("tanhLq.jl")

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
        stop::Bool,                     # "true" if routine stops at zero symdrome
        decayfactor::Float64,           # RBP decay factor
        listsizes::Vector{Int},         # sizes of the list for List-RBP algorithm
        rgn_seed::Int,                  # random seed to generate noise and message
        test::Bool,                     # if "true", perform test algorithm
        printtest::Bool;                # if "true", print test algorithm results
        msgtest=nothing,                # for testing specific msg
        noisetest=nothing               # for testing specific noise
    )::Tuple{Union{Matrix{Float64},Nothing},Union{Matrix{Float64},Nothing},Vector{Int},Vector{Int}}
    
################################## CONSTANTS ###################################

    # Parity-Check Matrix dimensions
    M,N = size(H)

    # 2*LS in NR5G (determines the number of initial punctured bits)
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
    Lf = (bptype != "MKAY") ? Vector{Float64}(undef,N) : Matrix{Float64}(undef,N,2)

    if twoLs > 0
        for i in 1:twoLs
            Lf[i] = 0.0
        end
    end

    # received signal
    signal = Vector{Float64}(undef,G)

    # bit-error
    biterror = Vector{Bool}(undef,N) 

    # Lq -> matrix of Nv messages (N x M)
    # Lr -> matrix of Nc messages (N x M)
    # if algorithm == "MKAY" the are different matrices for bit = 0 and bit = 1
    Lq = (bptype != "MKAY") ? Matrix{Float64}(undef,M,N) : Array{Float64,3}(undef,M,N,2)
    Lr = (bptype != "MKAY") ? Matrix{Float64}(undef,M,N) : Array{Float64,3}(undef,M,N,2)
    
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

    if algorithm == "RPD" || test
        syndrome = Vector{Bool}(undef,M)
    end

    if algorithm == "RPD"
        f = Vector{Int}(undef,N)
        oldLr = Matrix{Float64}(undef,M,N)
    end


############################### RBP PREALLOCATIONS #############################
    
    # RBP switchs 
    RBP_based = algorithm != "Flooding" && algorithm != "LBP" && algorithm != "RPD"

    if RBP_based

        # Consensus
        RBP_consensus = algorithm == "C-RBP"    || 
                        algorithm == "C&R-RBP"  ||
                        algorithm == "C&DR-RBP" 

        # Algorithms that used the RBP residual decay strategy
        RBP_decay = algorithm == "RD-RBP"   || 
                    algorithm == "List-RBP" || 
                    RBP_consensus

        # Algorithm that uses the same RBP routine ("RBP.jl")
        RBP_routine = algorithm == "RBP"      || 
                      algorithm == "RD-RBP"   || 
                      RBP_consensus

        # Algorithms that uses the default RBP initialization
        RBP_init = RBP_routine || algorithm == "NW-RBP" || algorithm == "SVNF" 

        # Variables common to all RBP based algorithms
        newLr = Matrix{Float64}(undef,M,N)
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

        if RBP_routine  
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
        elseif algorithm == "D-SVNF"
            u = zeros(Bool,N)    
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
            newLv = zeros(N)
            Lv = Vector{Float64}(undef,N)
            C = Vector{Bool}(undef,N)
            upc = zeros(Int,N)
        end
    end

################################## MAIN LOOP ###################################
    
    # for trial in 1:trials
    @fastmath @inbounds for trial in 1:trials

        ### 1) generate the random message
        generate_message!(msg,rgn,msgtest)

        ### 2) generate the cword
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
        received_signal!(signal,cword,G,twoLs,stdev,rgn,noisetest)

        # print info if in test algorithm
        if test && printtest
            println("Trial #$trial:")
            print_test("Message",msg)
            print_test("Transmitted codeword",cword[twoLs+1:end])
        end
        
        ### 4) reset simulation variables
        if test
            syndrome .= true
        end
        decoded .= false
        resetmatrix!(Lr,Nv,0.0)

        ### 5) init the LLR priors
        calc_Lf!(Lf,twoLs,signal,variance)
        if bptype == "TABL"
            # scale for table
            Lf[twoLs+1:end] .*= SIZE_PER_RANGE
        end
        # initial estimate of the bits (0-th iteration)
        for i in eachindex(bitvector)
            bitvector[i] = signbit(Lf[i])
        end

        ### 6) init the Lq matrix
        if algorithm == "RPD"
            init_Lq!(Lq,Lf,Nv,0.0)
        else
            init_Lq!(Lq,Lf,Nv,0.0)
        end
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
                init_residues!(Lq,Lr,Nc,phi,newLr,alpha,Residues,msum_factor)
            elseif algorithm == "List-RBP"
                list .= 0.0
                coords .= 0
                local_list .= 0.0
                local_coords .= 0
                init_list!(Lq,Lr,Nc,phi,msum_factor,newLr,Factors,inlist,Residues,list,
                                                                    coords,listsize)  
            elseif algorithm == "VC-RBP"
                for vj in eachindex(Nv)
                    alp = 0.0
                    for ci in Nv[vj]
                        li = LinearIndices(Lq)[ci,vj]
                        residue = abs(Lq[li])
                        if residue > alp
                            alp = residue
                        end
                        Residues[li] = residue
                    end
                    alpha[vj] = alp
                end
            elseif algorithm == "OV-RBP"
                for ci in eachindex(Nc)
                    Nci = Nc[ci]
                    for vj in Nci
                        newLr[ci,vj] = calc_Lr(Nci,ci,vj,Lq,msum_factor) 
                    end
                end
                size_C = 0
                C .= false
                for vj in eachindex(Nv)
                    Nvj = Nv[vj]
                    oldlv = Lf[vj]
                    newlv = calc_Ld(vj,Nvj,Lf,newLr)
                    newLv[vj] = newlv
                    Residues[vj] = abs(oldlv - newlv)
                    Lv[vj] = oldlv
                    if sign(oldlv)*sign(newlv) < 0
                        size_C += 1
                        C[vj] = true
                    end
                end

                max_upc = 0
                for vj in eachindex(Nv)
                    if C[vj]
                        count = 0
                        for ci in Nv[vj]
                            if _calc_syndrome(bitvector,Nc[ci])
                                count += 1
                            end
                        end
                        if count > max_upc
                            max_upc = count
                        end
                        upc[vj] = count
                    end
                end    
            else
                # D-SVNF                
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
                if RBP_routine
                    if switch_C_DR && iter ≥ 3
                        switch_R = true
                        switch_C_DR = false
                        num_reps -= N
                    end             
                    rbp_not_converged = RBP!(
                        bitvector,
                        Lq,
                        Lr,
                        Lf,
                        Nc,
                        Nv,
                        phi,
                        msum_factor,
                        msum2,
                        num_reps,
                        newLr,
                        Residues,
                        alpha,
                        RBP_decay,
                        decayfactor,
                        Factors,
                        rbp_not_converged,
                        RBP_consensus,
                        switch_R
                        )
                elseif algorithm == "NW-RBP"
                    rbp_not_converged = NW_RBP!(
                        bitvector,
                        Lq,
                        Lr,
                        Lf,
                        Nc,
                        Nv,
                        phi,
                        msum_factor,
                        msum2,
                        num_reps,
                        newLr,
                        alpha,
                        rbp_not_converged                    
                    )
                elseif algorithm == "SVNF"
                    rbp_not_converged = SVNF!(
                        bitvector,
                        Lq,
                        Lr,
                        Lf,
                        Nc,
                        Nv,
                        phi,
                        msum_factor,
                        msum2,
                        num_reps,
                        newLr,
                        Residues,
                        rbp_not_converged,
                        twoLs,
                        N
                    )
                elseif algorithm == "D-SVNF"
                    rbp_not_converged = D_SVNF!(
                        bitvector,
                        Lq,
                        Lr,
                        Lf,
                        Nc,
                        Nv,
                        phi,
                        msum_factor,
                        msum2,
                        num_reps,
                        newLr,
                        Residues,
                        rbp_not_converged,
                        N,
                        u
                    )
                elseif algorithm == "List-RBP"
                    rbp_not_converged = List_RBP!(
                        bitvector,
                        Lq,
                        Lr,
                        Lf,
                        Nc,
                        Nv,
                        phi,
                        msum_factor,
                        msum2,
                        num_reps,
                        newLr,
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
                        Lq,
                        Lr,
                        Lf,
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
                    rbp_not_converged, max_upc, size_C = OV_RBP!(
                        bitvector,
                        Lq,
                        Lr,
                        Lf,
                        Nc,
                        Nv,
                        phi,
                        msum_factor,
                        msum2,
                        num_reps,
                        newLr,
                        Residues,
                        rbp_not_converged,                    
                        Lv,
                        newLv,
                        C,
                        upc,
                        max_upc,
                        size_C
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
                    Lq,
                    Lr,
                    Lf,
                    Nc,
                    Nv,
                    phi,
                    msum_factor
                    )
            elseif algorithm == "LBP"
                LBP!(
                    bitvector,
                    Lq,
                    Lr,
                    Lf,
                    Nc,
                    Nv,
                    phi,
                    msum_factor
                    )
            elseif algorithm == "RPD"
                RPD!(
                    bitvector,
                    Lq,
                    Lr,
                    Lf,
                    Nc,
                    Nv,
                    phi,
                    msum_factor,
                    oldLr,
                    syndrome,
                    f
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
                retr[ci,vj] = log.(Lr[ci,vj,1]) - log.(Lr[ci,vj,2])
                retq[ci,vj] = log.(Lq[ci,vj,1]) - log.(Lq[ci,vj,2])
            end
        end
        Lr = retr
        Lq = retq
    end

    if test
        return Lr, Lq, sum_decoded, sum_ber
    else
        return nothing, nothing, sum_decoded, sum_ber
    end

end