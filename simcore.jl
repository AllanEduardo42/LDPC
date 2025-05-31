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
        A::Int,
        K::Int,
        R::Float64,
        G::Int,
        g_CRC::Vector{Bool},
        ebn0::Float64,
        H::Matrix{Bool},
        H1::Matrix{Bool},
        L::Matrix{Bool},
        U::Matrix{Bool},
        Nc::Vector{Vector{Int}},
        Nv::Vector{Vector{Int}},
        E_H::Matrix{Int},
        protocol::String,
        liftsize::Int,
        mode::String,
        bptype::String,
        trials::Int,
        maxiter::Int,
        stop::Bool,
        decayfactor::Float64,
        listsizes::Vector{Int},
        relative::Bool,
        rgn_seed_noise::Int,
        rgn_seed_msg::Int,
        test::Bool,
        printtest::Bool;
        msgtest=nothing,
        noisetest=nothing 
    )::Tuple{Matrix{Float64},Matrix{Float64},Vector{Int},Vector{Int}}
    
################################## CONSTANTS ###################################

    # Parity-Check Matrix dimensions
    M,N = size(H)

    # 2*liftsize in NR5G (determines the number of initial punctured bits)
    # = 0 otherwise
    twoZc = 0

    # protocol constants
    if protocol == "WiMAX" || protocol == "NR5G"
        E_M, E_N = size(E_H)
        E_K = E_N - E_M 
        if protocol == "WiMAX"            
            E_B = E_M  
            S = 0   
        elseif protocol == "NR5G"
            twoZc = 2*liftsize
            E_B = 4
            S = E_K*liftsize - K
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

################################# PREALLOCATIONS ###############################

    msg = Vector{Bool}(undef,A)  
    
    b = Vector{Bool}(undef,K)

    cword = Vector{Bool}(undef,N)

    if protocol == "PEG"
        Cw = cword
        w = Vector{Bool}(undef,M)
    else
        Cw = Matrix{Bool}(undef,liftsize,E_N) # info bits + CRC + filler bits + parity bits
        sw = Vector{Bool}(undef,liftsize)
        circ_aux = Vector{Bool}(undef,liftsize)    
        auxi = Vector{Bool}(undef,liftsize)   
        W = Matrix{Bool}(undef,liftsize,E_B)    
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

    if twoZc > 0
        if relative
            for i in 1:twoZc
                Lf[i] = 10/INFFLOAT
            end
        else
            for i in 1:twoZc
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

    # Set variables aux and signs depending on the BP type (for dispatching)
    if bptype == "FAST" || bptype == "MKAY"
        aux = Vector{Float64}(undef,N)
        signs = nothing
    elseif bptype == "TABL" 
        aux = Vector{Float64}(undef,N)
        signs = Vector{Bool}(undef,N)
    elseif bptype == "MSUM"
        aux = Vector{Float64}(undef,N)
        signs = Vector{Bool}(undef,N)
    else # bytype == TANH
        aux = nothing
        signs = nothing
    end

    phi = (bptype == "TABL") ? lookupTable() : nothing

    # RBP preallocations
    if mode == "RBP"
        newLr = Matrix{Float64}(undef,M,N)
        Factors = Matrix{Float64}(undef,M,N)
        resetmatrix!(Factors,Nv,1.0)
        residues = Vector{Float64}(undef,num_edges)
        coords = Matrix{Int}(undef,3,num_edges)
        rbpmatrix = Matrix{Int}(undef,M,N)
        e = 0
        for ci in eachindex(Nc)
            for vj in Nc[ci]
                e += 1
                coords[1,e] = ci
                coords[2,e] = vj
                li = LinearIndices(rbpmatrix)[ci,vj]
                coords[3,e] = li
                rbpmatrix[li] = e
            end
        end
        local_residues = nothing
        local_coords = nothing
    elseif mode == "List-RBP"
        newLr = Matrix{Float64}(undef,M,N)
        residues = Vector{Float64}(undef,listsizes[1]+1)
        Factors = Matrix{Float64}(undef,M,N)
        resetmatrix!(Factors,Nv,1.0)
        coords = Matrix{Int}(undef,3,listsizes[1]+1)
        rbpmatrix = Matrix{Bool}(undef,M,N)
        if listsizes[2] == 1
            local_residues = Vector{Float64}(undef,listsizes[1]+1)
            local_coords =  Matrix{Int}(undef,3,listsizes[1]+1)
        else
            local_residues = Vector{Float64}(undef,listsizes[2]+1)
            local_coords = Matrix{Int}(undef,3,listsizes[2]+1)
        end
    elseif mode == "NW-RBP"
        newLr = Matrix{Float64}(undef,M,N)
        residues = Vector{Float64}(undef,M)
        Factors = ones(M)   
    elseif mode == "VN-RBP" || mode == "VN-RBP2"
        newLr = Matrix{Float64}(undef,M,N)
        residues = Vector{Float64}(undef,N)
        Factors = ones(N)
        if mode == "VN-RBP2"
            mode2 = true
        else
            mode2 = false
        end
    elseif mode == "List-VN-RBP"
        newLr = Matrix{Float64}(undef,M,N)
        residues = Vector{Float64}(undef,listsizes[1]+1)
        Factors = ones(N)
        coords = Vector{Int}(undef,listsizes[1]+1)
        inlist = zeros(Bool,N)
        if listsizes[2] == 1
            local_residues = Vector{Float64}(undef,listsizes[1]+1)
            local_coords = Vector{Int}(undef,listsizes[1]+1)
        else
            local_residues = Vector{Float64}(undef,listsizes[2]+1)
            local_coords = Vector{Int}(undef,listsizes[2]+1)
        end
    end

################################## MAIN LOOP ###################################

    @inbounds for trial in 1:trials

        # 1) generate the random message
        generate_message!(msg,rgn_msg,msgtest)

        # 2) generate the cword
        append_CRC!(Cw,b,msg,g_CRC,A,K)
        if protocol == "PEG"
            encode_LDPC_LU!(Cw,H1,w,L,U,M,K)
            for i in 1:M
                Cw[K+i] = w[i]
            end
        else
            encode_LDPC_BG!(cword,Cw,W,sw,auxi,circ_aux,K,N,E_H,E_M,E_K,E_B,S,liftsize)
        end
        # verify the encoding
        if test
            _gf2_mat_mult!(syndrome,H,cword,M,N)
            if iszero(syndrome)
                println("Encoding error")
            end
        end

        # 3) sum the noise to the modulated cword to produce the received signal
        received_signal!(signal,cword,G,twoZc,stdev,rgn_noise,noisetest)

        # print info if in test mode
        if test && printtest
            println("Trial #$trial:")
            print_test("Message",msg)
            print_test("Transmitted cword",cword[twoZc+1:end])
        end
        
        # 4) reset simulation variables
        syndrome .= true
        decoded .= false
        resetmatrix!(Lr,Nv,0.0)
        if mode == "List-RBP" || mode == "List-VN-RBP"
            residues .= 0.0
            coords .= 0
            local_residues .= 0.0
            local_coords .= 0
            if mode == "List-RBP"
                resetmatrix!(rbpmatrix,Nv,false)
            end
        end

        # 5) init the LLR priors
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
        elseif mode == "VN-RBP" || mode == "VN-RBP2"
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
            elseif mode == "VN-RBP" || mode == "VN-RBP2"
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
                    rbp_not_converged,
                    mode2
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

            for vj in 1:N
                biterror[vj] = bitvector[vj] ≠ cword[vj]
            end
            
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

        # if all the residues are zero
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