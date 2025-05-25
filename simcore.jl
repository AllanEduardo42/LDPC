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
include("NW_RBP.jl")
include("./RBP functions/calc_all_residues.jl")
include("./RBP functions/calc_all_residues_NW.jl")
include("./RBP functions/calc_all_residues_VN.jl")
include("update_Lq.jl")
include("update_Lr.jl")
include("NR LDPC/NR_LDPC_functions.jl")
include("calc_parity_bits.jl")

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
        rgn_seed_msg::Integer;
        test=false,
        printtest=false,
        msgtest=nothing,
        noisetest=nothing        
    )

    if mode == "RBP" || mode == "List-RBP" || mode == "VN-RBP" || mode == "NW-RBP"
        RBP = true
    else
        RBP = false
    end
    
################################## CONSTANTS ###################################

    # Parity-Check Matrix dimensions
    M,N = size(H)

    twoZc = 0

    # constants for IEEE80216e
    if protocol == "PEG"
        Cw = Vector{Bool}(undef,N)
        z = Vector{Bool}(undef,M)
        v = Vector{Bool}(undef,M)
        w = Vector{Bool}(undef,M)
        cword = Cw
    else
        cword = Vector{Bool}(undef,N)
        E_M, E_N = size(E_H)
        E_K = E_N - E_M 
        Cw = Matrix{Bool}(undef,liftsize,E_N) # info bits + CRC + filler bits + parity bits
        Sw = Vector{Bool}(undef,liftsize)
        if protocol == "WiMAX"            
            E_B = E_M
            K = K_prime       
            P = 0   
        elseif protocol == "NR5G"
            E_B = 4
            twoZc = 2*liftsize
            K = E_K*liftsize
            S = K - K_prime             # quantity of filler bits
            P = liftsize*E_N - S - N    # quantity of removed parity bits (rate matching)
        end        
        W = Matrix{Bool}(undef,liftsize,E_B)
        Z = Matrix{Bool}(undef,liftsize,E_B)
    end
    # number of edges in the graph
    num_edges = sum(H)

    # transform EbN0 in standard deviations
    variance = exp10.(-ebn0/10) / (2*R)
    stdev = sqrt.(variance)

    # Set the random seeds
    rgn_noise = Xoshiro(rgn_seed_noise)
    rgn_msg = Xoshiro(rgn_seed_msg)

################################# PREALLOCATIONS ###############################

    msg = Vector{Bool}(undef,A)

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

    # prior llr (if mode == "MKAY" just the prior probabilities)
    if protocol == "NR5G" && relative
        a = 10/INFFLOAT
    else
        a = 0.0
    end
    Lf = (bptype != "MKAY") ? a*ones(N) : 0.5*ones(N,2)

    # noise
    noise = Vector{Float64}(undef,G)

    # received signal
    signal = Vector{Float64}(undef,G)

    # bit-error
    biterror = Vector{Bool}(undef,N) 

    # Lq -> matrix of Nv messages (N x M)
    # Lr -> matrix of Nc messages (M x N)
    # if mode == "MKAY" the are matrices for bit = 0 and bit = 1
    Lq = (bptype != "MKAY") ? zeros(M,N) : zeros(M,N,2)

    Lr = (bptype != "MKAY") ? zeros(M,N) : zeros(M,N,2)
    
    # Set variables aux and signs depending on the BP type (used for dispatch)
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

    phi = (bptype == "TABL") ? PHI : nothing

    # Set other variables that depend on the mode
    newLr = RBP ? H*0.0 : nothing
    Factors = RBP ? 1.0*H  : nothing
    
    # RBP modes
    if mode == "RBP"
        residues = zeros(num_edges)
        coords = zeros(Int,3,num_edges)
        localresidues = nothing
        localcoords = nothing
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
    elseif mode == "NW-RBP"
        residues = zeros(M)
        Factors = ones(M)
    elseif mode == "List-RBP"
        residues = zeros(listsizes[1]+1)
        coords = zeros(Int,3,listsizes[1]+1)
        rbpmatrix = Matrix(false*H)
        if listsizes[2] == 1
            localresidues = zeros(listsizes[1]+1)
            localcoords = zeros(Int,3,listsizes[1]+1)
        else
            localresidues = zeros(listsizes[2]+1)
            localcoords = zeros(Int,3,listsizes[2]+1)
        end
    elseif mode == "VN-RBP"
        residues = zeros(N)
        Factors = ones(N)
    end

################################## MAIN LOOP ###################################

    @inbounds for trial in 1:trials

        # 1) generate the random message
        generate_message!(msg,rgn_msg,msgtest)

        # 2) generate the cword
        Cw[1:A] = msg
        Cw[A+1:end] .= false
        _,Cw[A+1:K_prime] = divide_poly(Cw[1:K_prime],g_CRC)
        if protocol == "PEG"  
            z .= H[:,1:K_prime]*Cw[1:K_prime]
            gf2_solve_LU!(v,L,z)
            gf2_solve_LU!(w,U,v)
            Cw[K_prime+1:end] = w
        else
            calc_parity_bits!(Cw,W,Sw,Z,E_H,E_M,E_K,E_B)
            cword[1:K_prime] = Cw[1:K_prime]        # systematic bits
            cword[K_prime+1:end] = Cw[K+1:end-P]    # parity bits 
        end

        if !iszero(H*cword)
            println("ERROR")
        end

        # 3) Modulation of the cword
        @. signal = 2*cword[twoZc+1:end] - 1

        # 4) sum the noise to the modulated cword to produce the received signal
        received_signal!(signal,noise,stdev,rgn_noise,noisetest)

        # 5) print info in test mode
        if test && printtest
            println("Trial #$trial:")
            print_test("Message",msg)
            print_test("Transmitted cword",cword[twoZc+1:end])
        end
        
        # 6) reset simulation variables
        syndrome .= true
        decoded .= false
        resetmatrix!(Lr,Nv,0.0)
        if mode == "List-RBP"
            residues .= 0.0
            coords .= 0
            localresidues .= 0.0
            localcoords .= 0
            resetmatrix!(rbpmatrix,Nv,false)
        end

        # 7) init the LLR priors
        calc_Lf!(Lf,twoZc,signal,variance)
        if bptype == "TABL"
            # scale for table
            Lf[twoZc+1:end] .*= SIZE_PER_RANGE
        end
        for i in eachindex(bitvector)
            bitvector[i] = signbit(Lf[i])
        end
        
        # 0-th iter
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
        end
        
        # 10) BP routine
        bp_not_converged = true

        iter = 0
        while iter < maxiter && bp_not_converged
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
                bp_not_converged = RBP!(
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
                    localresidues,
                    localcoords,
                    listsizes,
                    relative,
                    bp_not_converged
                    )
                # reset factors
                resetmatrix!(Factors,Nv,1.0)
            elseif mode == "NW-RBP"
                bp_not_converged = NW_RBP!(
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
                    bp_not_converged
                )
                # reset factors
                Factors .= 1.0
            elseif mode == "VN-RBP"
                bp_not_converged = VN_RBP!(
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
                    bp_not_converged
                    )
                # reset factors
                Factors .= 1.0
            end            
    
            calc_syndrome!(syndrome,bitvector,Nc)

            biterror .= (bitvector .≠ cword)
    
            if test
                if printtest
                    print_test("Bit error",biterror)   
                    println("Bit error rate: $(sum(biterror))/$N")
                    print_test("Syndrome",syndrome)  
                    println("Syndrome rate: $(sum(syndrome))/$M")
                    println() 
                    if !bp_not_converged
                        println("#### BP has converged at iteration $iter ####")
                    end
                end
                if iszero(syndrome)
                    if iszero(biterror)
                        display("everyting is fine")
                    else
                        display("wrong decoding")
                        break
                    end
                end
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

        if !bp_not_converged
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
            return retr, retq
        else
            return Lr, Lq
        end
    else
        return sum_decoded, sum_ber
    end

end