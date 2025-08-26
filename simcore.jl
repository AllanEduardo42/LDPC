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
include("VN_LBP.jl")
include("RBP.jl")
include("E_NV_RBP.jl")
include("List_RBP.jl")
include("List_VN_RBP.jl")
include("SVNF.jl")
include("VN_RBP.jl")
include("LD_RBP.jl")
include("NW_RBP.jl")
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
        mode::String,                   # BP algorithm (flooding, RBP etc.)
        bptype::String,                 # Type of BP implementation (FAST, TANH etc.)
        trials::Int,                    # Number of trials
        maxiter::Int,                   # Maximum number of BP iterations
        stop::Bool,                     # "true" if routine stops at zero symdrome
        decayfactor::Float64,                     # RBP decay factor
        listsizes::Vector{Int},         # sizes of the list for List-RBP mode
        rgn_seed::Int,                  # random seed to generate noise and message
        test::Bool,                     # if "true", perform test mode
        printtest::Bool;                # if "true", print test mode results
        msgtest=nothing,                # for testing specific msg
        noisetest=nothing               # for testing specific noise
    )::Tuple{Matrix{Float64},Matrix{Float64},Vector{Int},Vector{Int}}
    
################################## CONSTANTS ###################################

    # Parity-Check Matrix dimensions
    M,N = size(H)

    # 2*LS in NR5G (determines the number of initial punctured bits)
    # = 0 otherwise
    twoLs = 0

    if mode == "RBP relative"
        relative = true
        mode = "RBP"
    else
        relative = false
    end

    # protocol constants
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

    # syndrome
    syndrome = Vector{Bool}(undef,M)

    # prior llr (if mode == "MKAY" it is just the prior probabilities)
    Lf = (bptype != "MKAY") ? Vector{Float64}(undef,N) : Matrix{Float64}(undef,N,2)

    if twoLs > 0
        if relative
            for i in 1:twoLs
                Lf[i] = 0.001
            end
        else
            for i in 1:twoLs
                Lf[i] = 0.0
            end
        end
    end

    # received signal
    signal = Vector{Float64}(undef,G)

    # bit-error
    biterror = Vector{Bool}(undef,N) 

    # Lq -> matrix of Nv messages (N x M)
    # Lr -> matrix of Nc messages (N x M)
    # if mode == "MKAY" the are different matrices for bit = 0 and bit = 1
    Lq = (bptype != "MKAY") ? Matrix{Float64}(undef,M,N) : Array{Float64,3}(undef,M,N,2)
    Lr = (bptype != "MKAY") ? Matrix{Float64}(undef,M,N) : Array{Float64,3}(undef,M,N,2)

    # Set variables signs depending on the BP type (for dispatching)
    if bptype == "TABL" || bptype == "MSUM"
        signs = Vector{Bool}(undef,N)
    else
        signs = nothing
    end

    phi = (bptype == "TABL") ? lookupTable() : nothing

############################### RBP PREALLOCATIONS #############################
    
    LIST = mode == "List-RBP" || mode == "List-VN-RBP"
    RBP_mode = mode[end-2:end] == "RBP" || mode == "SVNF"
    if RBP_mode
        newLr = Matrix{Float64}(undef,M,N)
        Residues = Matrix{Float64}(undef,M,N)
        NOT_SVNF = mode ≠ "SVNF"
        NOT_NW = mode ≠ "NW-RBP"
        if NOT_SVNF
            alpha = Vector{Float64}(undef,M)
        end
        if NOT_SVNF && NOT_NW
            Factors = Matrix{Float64}(undef,M,N)
            resetmatrix!(Factors,Nv,1.0)
        end 
        switch_C_VN = false  
        if mode == "TW-RBP" 
            num_reps = num_edges - N
            two_ways = true
            consensus = false
            C_VN = false      
        elseif mode == "C-RBP"
            num_reps = num_edges
            two_ways = false
            consensus = true
            C_VN = false                     
        elseif mode == "C-VN-RBP"  
            num_reps = num_edges
            two_ways = false
            consensus = true
            C_VN = true
        elseif mode == "RBP"
            num_reps = num_edges
            two_ways = false
            consensus = false
            C_VN = false
        end
        if LIST
            listsize = listsizes[1]
            list = Vector{Float64}(undef,listsize+1)
            list[end] = 0.0
            coords = Matrix{Int}(undef,3,listsize+1)
            inlist = zeros(Bool,M,N)
            if mode == "List-RBP"
                listsize2 = listsizes[2]
                local_list = Vector{Float64}(undef,listsize2+1)
                local_coords = Matrix{Int}(undef,3,listsize2+1)
            end            
        end       
    # elseif mode == "LD-RBP"
    #     newLr = Matrix{Float64}(undef,M,N)
    #     alpha = Vector{Float64}(undef,N)
    #     Factors = ones(N)
    end

    

################################## MAIN LOOP ###################################

    @inbounds for trial in 1:trials

        if mode == "List-VN-RBP"
            listsize = 2
        end

        # 1) generate the random message
        generate_message!(msg,rgn,msgtest)

        # 2) generate the cword
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
        # verify the encoding
        if test
            _gf2_mat_mult!(syndrome,H,cword,M,N)
            if !iszero(syndrome)
                throw(error("encoding error"))
            end
        end

        # 3) sum the noise to the modulated cword to produce the received signal
        received_signal!(signal,cword,G,twoLs,stdev,rgn,noisetest)

        # print info if in test mode
        if test && printtest
            println("Trial #$trial:")
            print_test("Message",msg)
            print_test("Transmitted codeword",cword[twoLs+1:end])
        end
        
        # 4) reset simulation variables
        syndrome .= true
        decoded .= false
        resetmatrix!(Lr,Nv,0.0)

        if LIST
            list .= 0.0
            coords .= 0
            resetmatrix!(inlist,Nv,false)
            if mode == "List-RBP"
                local_list .= 0.0
                local_coords .= 0
            end
        end

        # 5) init the LLR priors
        calc_Lf!(Lf,twoLs,signal,variance)
        if bptype == "TABL"
            # scale for table
            Lf[twoLs+1:end] .*= SIZE_PER_RANGE
        end
        # initial estimate of the bits (0-th iteration)
        for i in eachindex(bitvector)
            bitvector[i] = signbit(Lf[i])
        end
        
        # print info if in test mode
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
        
        # 8) init the Lq matrix
        init_Lq!(Lq,Lf,Nv,signs)

        # 9) init the RBP methods
        if mode == "RBP" || mode == "TW-RBP" || mode == "C-RBP" || mode == "VN-RBP" || mode == "C-VN-RBP"
            init_RBP!(Lq,Lr,Nc,signs,phi,newLr,alpha,Residues)
        elseif LIST
            init_list!(Lq,Lr,Nc,signs,phi,newLr,Factors,inlist,Residues,list,
                                                                coords,listsize)
        elseif mode == "SVNF"
            init_SVNF!(Lq,Lr,Nc,signs,phi,newLr,Residues)        
        elseif mode == "NW-RBP"
            init_NW_RBP!(Lq,Lr,Nc,signs,phi,newLr,alpha)
        elseif mode == "LD-RBP"
            init_LD_RBP!(Lq,Nc,signs,phi,newLr,alpha,Nv)
        end
  
        # 10) BP routine
        rbp_not_converged = true

        iter = 0
        while iter < maxiter && rbp_not_converged
            iter += 1

            if test && printtest  
                println("### Iteration #$iter ###")
            end
    
            if mode == "Flooding"
                flooding!(
                    bitvector,
                    Lq,
                    Lr,
                    Lf,
                    Nc,
                    Nv,
                    signs,
                    phi
                    )
            elseif mode == "LBP"
                LBP!(
                    bitvector,
                    Lq,
                    Lr,
                    Lf,
                    Nc,
                    Nv,
                    signs,
                    phi
                    )
            elseif mode == "VN-LBP"
                VN_LBP!(
                    bitvector,
                    Lq,
                    Lr,
                    Lf,
                    Nc,
                    Nv)
            elseif mode == "RBP" || mode == "TW-RBP" || mode == "C-RBP" || mode == "C-VN-RBP"
                if C_VN && iter ≥ 4
                    switch_C_VN = true
                    C_VN = false
                    num_reps -= N
                end
                rbp_not_converged = RBP!(
                    bitvector,
                    Lq,
                    Lr,
                    Lf,
                    Nc,
                    Nv,
                    signs,
                    phi,
                    decayfactor,
                    num_reps,
                    newLr,
                    alpha,
                    Residues,
                    Factors,
                    rbp_not_converged,
                    two_ways,
                    consensus,
                    switch_C_VN
                    )
                # reset factors
                resetmatrix!(Factors,Nv,1.0)
            elseif mode == "VN-RBP"
                rbp_not_converged = VN_RBP!(
                    bitvector,
                    Lq,
                    Lr,
                    Lf,
                    Nc,
                    Nv,
                    signs,
                    phi,
                    decayfactor,
                    num_edges-N,
                    newLr,
                    alpha,
                    Residues,
                    Factors,
                    rbp_not_converged
                    )
                # reset factors
                resetmatrix!(Factors,Nv,1.0)
            elseif mode == "List-RBP"
                rbp_not_converged = List_RBP!(
                    bitvector,
                    Lq,
                    Lr,
                    Lf,
                    Nc,
                    Nv,
                    signs,
                    phi,
                    decayfactor,
                    num_edges,
                    newLr,
                    Factors,
                    coords,
                    inlist,
                    Residues,
                    list,
                    local_list,
                    local_coords,
                    listsize,
                    listsize2,
                    rbp_not_converged
                )
                # reset factors
                resetmatrix!(Factors,Nv,1.0)
            elseif mode == "SVNF"
                rbp_not_converged = SVNF!(
                    bitvector,
                    Lq,
                    Lr,
                    Lf,
                    Nc,
                    Nv,
                    signs,
                    phi,
                    N,
                    num_edges,
                    newLr,
                    Residues,
                    rbp_not_converged,
                    twoLs
                )
            elseif mode == "NW-RBP"
                rbp_not_converged = NW_RBP!(
                    bitvector,
                    Lq,
                    Lr,
                    Lf,
                    Nc,
                    Nv,
                    signs,
                    phi,
                    # decayfactor,
                    M,
                    newLr,
                    # Factors,
                    alpha,
                    rbp_not_converged
                )
                # reset factors
                # Factors .= 1.0
            elseif mode == "List-VN-RBP"
                if listsize < listsizes[1]
                    listsize = 2^(cld(iter,2)+1)
                    # listsize2 = listsize
                end
                rbp_not_converged = List_VN_RBP!(
                    bitvector,
                    Lq,
                    Lr,
                    Lf,
                    Nc,
                    Nv,
                    signs,
                    phi,
                    decayfactor,
                    num_edges-N,
                    newLr,
                    Residues,
                    list,
                    Factors,
                    coords,
                    inlist,
                    listsize,
                    rbp_not_converged
                    )
                # reset factors
                resetmatrix!(Factors,Nv,1.0)
            elseif mode == "LD-RBP"
               rbp_not_converged = LD_RBP!(
                    bitvector,
                    Lq,
                    Lr,
                    Lf,
                    Nc,
                    Nv,
                    signs,
                    phi,
                    decayfactor,
                    num_edges-N,
                    newLr,
                    Factors,
                    alpha,
                    rbp_not_converged
                    )
                # reset factors
                Factors .= 1.0
            end            

            for vj in 1:N
                biterror[vj] = bitvector[vj] ≠ cword[vj]
            end
            
            # print info if in test mode
            if test
                calc_syndrome!(syndrome,bitvector,Nc)
                if printtest
                    print_test("Syndrome",syndrome)  
                    println("Syndrome rate: $(sum(syndrome))/$M")
                    print_test("Bit error",biterror)   
                    println("Bit error rate: $(sum(biterror))/$N")
                    println()
                    if iszero(biterror)                        
                        decoded[iter] = true
                        rbp_not_converged = false
                    end
                    if !rbp_not_converged
                        println("#### BP has converged at iteration $iter ####")
                    end
                end
            else
                zero_syn = iszerosyndrome(bitvector,Nc)
                if zero_syn
                    if iszero(biterror)                        
                        decoded[iter] = true
                        rbp_not_converged = false
                    end
                    if stop
                        if iter < maxiter
                            for i = iter+1:maxiter
                                decoded[i] = decoded[iter]
                            end
                        end
                        break
                    end
                end               
                ber[iter] = sum(biterror) 
            end
        end

        # if all the Residues are zero
        if !rbp_not_converged
            for i = iter+1:maxiter
                decoded[i] = decoded[iter]
                ber[i] = ber[iter]
            end
        end

        # bit error rate
        for i=1:maxiter
            sum_ber[i] += ber[i]
            sum_decoded[i] += decoded[i]
        end

        if mode == "C-VN-RBP"
            switch_C_VN = false
            C_VN = true
            num_reps = num_edges
        end
    end

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

    return Lr, Lq, sum_decoded, sum_ber

end