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

function
    simcore(
        A::Integer,
        R::Rational,
        G::Integer,
        ebn0::AbstractFloat,
        H::Matrix{Bool},
        L::Union{Nothing,Matrix{Bool}},
        U::Union{Nothing,Matrix{Bool}},
        Nc::Vector{Vector{T}} where {T<:Integer},
        Nv::Vector{Vector{T}} where {T<:Integer},
        E_H::Union{Nothing,Matrix{<:Integer}},
        protocol::String,
        Zf::Integer,
        Nr_ldpc_data::Union{Nothing,nr_ldpc_data},
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
    if protocol == "WiMAX"
        E_M, E_N = size(E_H)
        E_K = E_N - E_M
        Cw = zeros(Bool,Zf,E_N+1)
        W = zeros(Bool,Zf,E_M)
        Sc = zeros(Bool,Zf)
    elseif protocol == "NR5G"
        B = Nr_ldpc_data.B
        K_prime = Nr_ldpc_data.K_prime
        K = Nr_ldpc_data.K
        Zc = Nr_ldpc_data.Zc
        bg = Nr_ldpc_data.bg
        N_cb = Nr_ldpc_data.N_cb
        E_r = Nr_ldpc_data.E_r
        k0 = Nr_ldpc_data.k0
        g_CRC = Nr_ldpc_data.g_CRC
        P_Zc = Nr_ldpc_data.P_Zc
        
        twoZc = 2*Zc
        W = Matrix{Bool}(undef,Zc,4)
        Z = Matrix{Bool}(undef,Zc,4)
        Sc = Vector{Bool}(undef,Zc)
        if bg == "1"
            Cw = Matrix{Bool}(undef,Zc,68-P_Zc)
            J = 22
            I = 46 - Int(P_Zc)
        else
            Cw = Matrix{Bool}(undef,Zc,52-P_Zc)
            J = 10
            I = 42 - Int(P_Zc)
        end
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

    cword = Vector{Bool}(undef,G)

    complete_cword = Vector{Bool}(undef,N)

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
    signal = zeros(G)

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
        if protocol == "NR5G"
            Cw[1:A] = msg
            Cw[A+1:B] .= false
            _,Cw[A+1:K_prime] = divide_poly(Cw[1:B],g_CRC)
            Cw[K_prime+1:end] .= false
            W .= false
            Sc .= false
            Z .= false
            NR_LDPC_parity_bits!(Cw,W,Sc,Z,E_H,I,J)
            rate_matching!(cword,Cw,twoZc,N_cb,E_r[1],k0,K,K_prime)
        elseif protocol == "PEG"
            LU_parity_bits!(cword,H,L,U,msg)
        elseif protocol == "WiMAX"
            Cw .= false
            W .= false
            Sc .= false
            Cw[1:Zf*E_K] = msg
            IEEE80216e_parity_bits!(Cw,W,Sc,Zf,E_H,E_M,E_K)
            cword .= Cw[1:end-Zf]
        end

        # 3) Modulation of the cword
        @. signal = 2*cword - 1

        # 4) sum the noise to the modulated cword to produce the received signal
        received_signal!(signal,noise,stdev,rgn_noise,noisetest)

        # 5) print info in test mode
        if test && printtest
            println("Trial #$trial:")
            print_test("msg",msg)
            print_test("cword",cword)
        end
        
        # 6) Include the punctured bits in the cword for biterror calculation
        if twoZc > 0
            complete_cword[1:twoZc] = msg[1:twoZc]
            complete_cword[twoZc+1:end] = cword
        else
            complete_cword .= cword
        end

        # 7) reset simulation variables
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

        # 8) init the LLR priors
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
            biterror .= (bitvector .≠ complete_cword)
            print_test("Bit error",biterror)   
            println("Bit error rate: $(sum(biterror))/$N")
            print_test("Syndrome",syndrome)  
            println("Syndrome rate: $(sum(syndrome))/$M")
            println()
        end
        
        # 9) init the Lq matrix
        init_Lq!(Lq,Lf,Nv)

        # 10) precalculate the residues for RBP
        if mode == "RBP" || mode == "List-RBP"
            calc_all_residues!(Lq,Lr,Nc,aux,signs,phi,newLr,Factors,rbpmatrix,
                                        residues,coords,listsizes,relative)
        elseif mode == "NW-RBP"
            calc_all_residues_NW!(Lq,Nc,aux,signs,phi,newLr,residues)
        elseif mode == "VN-RBP"
            calc_all_residues_VN!(Lq,Nc,aux,signs,phi,Lr,newLr,residues,Nv)
        end
        
        # BP routine
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

            biterror .= (bitvector .≠ complete_cword)
    
            if test && printtest
                print_test("Bit error",biterror)   
                println("Bit error rate: $(sum(biterror))/$N")
                print_test("Syndrome",syndrome)  
                println("Syndrome rate: $(sum(syndrome))/$M")
                println() 
                if !bp_not_converged
                    println("#### BP has converged at iteration $iter ####")
                end        
            else
                if iszero(syndrome)
                    if bitvector == complete_cword
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