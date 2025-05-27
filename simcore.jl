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
include("VN_RBP.jl")
include("List_VN_RBP.jl")
include("NW_RBP.jl")
include("./RBP functions/calc_all_residues.jl")
include("./RBP functions/calc_all_residues_NW.jl")
include("./RBP functions/calc_all_residues_VN.jl")
include("./RBP functions/calc_all_residues_list_VN.jl")
include("update_Lq.jl")
include("update_Lr.jl")
include("NR LDPC/NR_LDPC_functions.jl")
include("encode_LDPC.jl")

function
    simcore(
        A::Integer,
        K_prime::Integer,
        R::AbstractFloat,
        G::Integer,
        g_CRC::Vector{Bool},
        ebn0::AbstractFloat,
        H::Matrix{Bool},
        L::Union{Nothing,Matrix{Bool}},
        U::Union{Nothing,Matrix{Bool}},
        Nc::Vector{Vector{T}} where {T<:Integer},
        Nv::Vector{Vector{T}} where {T<:Integer},
        E_H::Union{Nothing,Matrix{<:Integer}},
        protocol::String,
        liftsize::Integer,
        mode::String,
        bptype::String,
        trials::Integer,
        maxiter::Integer,
        stop::Bool,
        decayfactor::AbstractFloat,
        listsizes::Vector{<:Integer},
        relative::Bool,
        rgn_seed_noise::Integer,
        rgn_seed_msg::Integer,
        test::Bool,
        printtest::Bool,
        msgtest::Union{Nothing,Vector{Bool}},
        noisetest::Union{Nothing,Vector{Bool}}   
    )
    
################################## CONSTANTS ###################################

    # Parity-Check Matrix dimensions
    M,N = size(H)

    # 2*liftsize in NR5G (determines the number of initial punctured bits)
    # =0 otherwise
    twoZc = 0

    # protocol constants
    if protocol == "WiMAX" || protocol == "NR5G"
        E_M, E_N = size(E_H)
        E_K = E_N - E_M 
        if protocol == "WiMAX"            
            E_B = E_M
            K = K_prime       
            P = 0   
        elseif protocol == "NR5G"
            twoZc = 2*liftsize
            E_B = 4
            K = E_K*liftsize
            P = liftsize*E_N - (K - K_prime) - N # quantity of removed parity bits (rate matching)
        end
    end

    # number of edges in the graph
    num_edges = sum(H)

    # transform EbN0 in standard deviations
    variance = exp10.(-ebn0/10) / (2*(R+16/G))
    stdev = sqrt.(variance)

    # Set the random seeds
    rgn_noise = Xoshiro(rgn_seed_noise)
    rgn_msg = Xoshiro(rgn_seed_msg)

    # if NR5G and relative, we must init the llr with not zero values due to puncturing
    if protocol == "NR5G" && relative
        init_Lf = 10/INFFLOAT
    else
        init_Lf = 0.0
    end

################################# PREALLOCATIONS ###############################

    msg = Vector{Bool}(undef,A)    

    cword = Vector{Bool}(undef,N)

    if protocol == "PEG"
        Cw = cword
        z = Vector{Bool}(undef,M)
        v = Vector{Bool}(undef,M)
        w = Vector{Bool}(undef,M)
    else
        Cw = Matrix{Bool}(undef,liftsize,E_N) # info bits + CRC + filler bits + parity bits
        sw = Vector{Bool}(undef,liftsize)       
        W = Matrix{Bool}(undef,liftsize,E_B)
        Z = Matrix{Bool}(undef,liftsize,E_B)
    end

    # frame error rate
    sum_decoded = zeros(Int,maxiter)
    decoded = Vector{Bool}(undef,maxiter)

    # bit error rate
    sum_ber = zeros(Int,maxiter)
    ber = zeros(Int,maxiter)

    # estimate
    bitvector = Vector{Bool}(undef,N)

    # syndrome
    syndrome = Vector{Bool}(undef,M)

    # prior llr (if mode == "MKAY" it is just the prior probabilities)
    Lf = (bptype != "MKAY") ? init_Lf*ones(N) : 0.5*ones(N,2)

    # noise
    noise = Vector{Float64}(undef,G)

    # received signal
    signal = Vector{Float64}(undef,G)

    # bit-error
    biterror = Vector{Bool}(undef,N) 

    # Lq -> matrix of Nv messages (N x M)
    # Lr -> matrix of Nc messages (N x M)
    # if mode == "MKAY" the are different matrices for bit = 0 and bit = 1
    Lq = (bptype != "MKAY") ? zeros(M,N) : zeros(M,N,2)
    Lr = (bptype != "MKAY") ? zeros(M,N) : zeros(M,N,2)
    
    # Set variables aux and signs depending on the BP type (for dispatching)
    if bptype == "FAST" || bptype == "MKAY"
        aux = zeros(N)
        signs = nothing
    elseif bptype == "TABL" 
        aux = zeros(N)
        signs = zeros(Bool,N)
    elseif bptype == "MSUM"
        aux = nothing
        signs = zeros(Bool,N)
    else # bytype == TANH
        aux = nothing
        signs = nothing
    end

    phi = (bptype == "TABL") ? lookupTable() : nothing

    # RBP preallocations
    if mode == "RBP"
        newLr = zeros(M,N)
        residues = zeros(num_edges)
        Factors = Float64.(H)
        coords = zeros(Int,3,num_edges)
        rbpmatrix = 0*H
        e = 0
        for m in axes(H,1)
            for n in Nc[m]
                e += 1
                coords[1,e] = m
                coords[2,e] = n
                li = LinearIndices(rbpmatrix)[m,n]
                coords[3,e] = li
                rbpmatrix[li] = e
            end
        end
        local_residues = nothing
        local_coords = nothing
    elseif mode == "List-RBP"
        newLr = zeros(M,N)
        residues = zeros(listsizes[1]+1)
        Factors = Float64.(H)
        coords = zeros(Int,3,listsizes[1]+1)
        rbpmatrix = Matrix(false*H)
        if listsizes[2] == 1
            local_residues = zeros(listsizes[1]+1)
            local_coords = zeros(Int,3,listsizes[1]+1)
        else
            local_residues = zeros(listsizes[2]+1)
            local_coords = zeros(Int,3,listsizes[2]+1)
        end
    elseif mode == "NW-RBP"
        newLr = zeros(M,N)
        residues = zeros(M)
        Factors = ones(M)       
    elseif mode == "VN-RBP"
        newLr = zeros(M,N)
        residues = zeros(N)
        Factors = ones(N)
    elseif mode == "List-VN-RBP"
        newLr = zeros(M,N)
        residues = zeros(listsizes[1]+1)
        Factors = ones(N)
        coords = zeros(Int,listsizes[1]+1)
        inlist = zeros(Bool,N)
        if listsizes[2] == 1
            local_residues = zeros(listsizes[1]+1)
            local_coords = zeros(Int,listsizes[1]+1)
        else
            local_residues = zeros(listsizes[2]+1)
            local_coords = zeros(Int,listsizes[2]+1)
        end
    end

################################## MAIN LOOP ###################################

    @inbounds for trial in 1:trials

        # 1) generate the random message
        generate_message!(msg,rgn_msg,msgtest)

        # 2) generate the cword
        append_CRC!(Cw,msg,g_CRC,A,K_prime)
        if protocol == "PEG"
            encode_LDPC_LU!(Cw,H,z,w,v,L,U,K_prime)
        else
            encode_LDPC_BG!(cword,Cw,Z,W,sw,K_prime,K,E_H,E_M,E_K,E_B,P)
        end        

        # verify the encoding
        if !iszero(H*cword)
            println("Encoding error")
        end

        # 3) Modulation of the cword in BPSK
        @. signal = 2*cword[twoZc+1:end] - 1

        # 4) sum the noise to the modulated cword to produce the received signal
        received_signal!(signal,noise,stdev,rgn_noise,noisetest)

        # print info if in test mode
        if test && printtest
            println("Trial #$trial:")
            print_test("Message",msg)
            print_test("Transmitted cword",cword[twoZc+1:end])
        end
        
        # 5) reset simulation variables
        syndrome .= true
        decoded .= false
        resetmatrix!(Lr,Nv,0.0)
        if mode == "List-RBP"
            residues .= 0.0
            coords .= 0
            local_residues .= 0.0
            local_coords .= 0
            resetmatrix!(rbpmatrix,Nv,false)
        end

        # 6) init the LLR priors
        calc_Lf!(Lf,twoZc,signal,variance)
        if bptype == "TABL"
            # scale for table
            Lf[twoZc+1:end] .*= SIZE_PER_RANGE
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
            print_test("Bit error",biterror)   
            println("Bit error rate: $(sum(biterror))/$N")
            print_test("Syndrome",syndrome)  
            println("Syndrome rate: $(sum(syndrome))/$M")
            println()
        end
        
        # 8) init the Lq matrix
        init_Lq!(Lq,Lf,Nv)

        # 9) precalculate the residues for RBP
        if mode == "RBP" || mode == "List-RBP"
            calc_all_residues!(Lq,Lr,Nc,aux,signs,phi,newLr,Factors,rbpmatrix,
                                        residues,coords,listsizes,relative)
        elseif mode == "NW-RBP"
            calc_all_residues_NW!(Lq,Nc,aux,signs,phi,newLr,residues)
        elseif mode == "VN-RBP"
            calc_all_residues_VN!(Lq,Nc,aux,signs,phi,Lr,newLr,residues,Nv)
        elseif mode == "List-VN-RBP"
            calc_all_residues_list_VN!(Lq,Nc,aux,signs,phi,Lr,newLr,residues,Nv,inlist,listsizes,coords)
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
                    aux,
                    signs,
                    phi)
            elseif mode == "LBP"
                LBP!(
                    bitvector,
                    Lq,
                    Lr,
                    Lf,
                    Nc,
                    Nv,
                    aux,
                    signs,
                    phi)
            elseif (mode == "RBP" || mode == "List-RBP")
                rbp_not_converged = RBP!(
                    bitvector,
                    Lq,
                    Lr,
                    Lf,
                    Nc,
                    Nv,
                    aux,
                    signs,
                    phi,
                    decayfactor,
                    num_edges,
                    newLr,
                    Factors,
                    coords,
                    rbpmatrix,
                    residues,
                    local_residues,
                    local_coords,
                    listsizes,
                    relative,
                    rbp_not_converged
                    )
                # reset factors
                resetmatrix!(Factors,Nv,1.0)
            elseif mode == "NW-RBP"
                rbp_not_converged = NW_RBP!(
                    bitvector,
                    Lq,
                    Lr,
                    Lf,
                    Nc,
                    Nv,
                    aux,
                    signs,
                    phi,
                    decayfactor,
                    M,
                    newLr,
                    Factors,
                    residues,
                    rbp_not_converged
                )
                # reset factors
                Factors .= 1.0
            elseif mode == "VN-RBP"
                rbp_not_converged = VN_RBP!(
                    bitvector,
                    Lq,
                    Lr,
                    Lf,
                    Nc,
                    Nv,
                    aux,
                    signs,
                    phi,
                    decayfactor,
                    num_edges,
                    newLr,
                    Factors,
                    residues,
                    rbp_not_converged
                    )
                # reset factors
                Factors .= 1.0
            elseif mode == "List-VN-RBP"
                rbp_not_converged = List_VN_RBP!(
                    bitvector,
                    Lq,
                    Lr,
                    Lf,
                    Nc,
                    Nv,
                    aux,
                    signs,
                    phi,
                    decayfactor,
                    num_edges,
                    newLr,
                    Factors,
                    residues,
                    rbp_not_converged,
                    inlist,
                    local_residues,
                    local_coords,
                    listsizes,
                    coords
                )
                # reset factors
                Factors .= 1.0
            end            
    
            calc_syndrome!(syndrome,bitvector,Nc)

            biterror .= (bitvector .≠ cword)
            
            # print info if in test mode
            if test
                if printtest
                    print_test("Bit error",biterror)   
                    println("Bit error rate: $(sum(biterror))/$N")
                    print_test("Syndrome",syndrome)  
                    println("Syndrome rate: $(sum(syndrome))/$M")
                    println() 
                    if !rbp_not_converged
                        println("#### BP has converged at iteration $iter ####")
                    end
                end
                # if iszero(syndrome)
                #     if iszero(biterror)
                #         display("everyting is fine")
                #     else
                #         display("wrong decoding")
                #         break
                #     end
                # end
            else
                if iszero(syndrome)
                    if iszero(biterror)
                        decoded[iter] = true
                    end
                    if stop
                        if iter < maxiter
                            decoded[iter+1:end] .= decoded[iter]
                        end
                        break
                    end
                end               
                ber[iter] = sum(biterror) 
            end
        end

        # if all the residues are zero
        if !rbp_not_converged
            decoded[iter+1:end] .= decoded[iter]
            ber[iter+1:end] .= ber[iter]
        end

        # bit error rate
        @. sum_ber += ber
        @. sum_decoded += decoded

    end

    if test
        if bptype == "MKAY"
            retr = zeros(M,N)
            retq = zeros(M,N)
            for m in eachindex(Nc)
                for n in Nc[m]
                    retr[m,n] = log.(Lr[m,n,1]) - log.(Lr[m,n,2])
                    retq[m,n] = log.(Lq[m,n,1]) - log.(Lq[m,n,2])
                end
            end
            Lr = retr
            Lq = retq
        end
        return Lr, Lq, nothing, nothing
    else
        return nothing, nothing, sum_decoded, sum_ber
    end
end