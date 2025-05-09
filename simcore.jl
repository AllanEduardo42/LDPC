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
include("Genius-RBP.jl")
include("VN_RBP.jl")
include("NW_RBP.jl")
include("./RBP functions/calc_all_residues.jl")
include("update_Lq.jl")
include("update_Lr.jl")

function
    simcore(
        A::Integer,
        R::AbstractFloat,
        ebn0::AbstractFloat,
        H::Matrix{Bool},
        G::Union{Nothing,Matrix{Bool}},
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

    if mode == "RBP" || mode == "List-RBP" || mode == "VN-RBP" || 
       mode == "Genius-RBP" || mode == "NW-RBP"
        RBP = true
    else
        RBP = false
    end
    
################################## CONSTANTS ###################################

    # Parity-Check Matrix dimensions
    M,N = size(H)
    
    # constant Zc of NR/5G
    if protocol == "NR5G"
        twoZc = 2*Nr_ldpc_data.Zc
    else
        twoZc = 0
    end

    # constants for IEEE80216e
    if protocol == "WiMAX"
        E_M, E_N = size(E_H)
        E_K = E_N - E_M
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
    Lf = (bptype != "MKAY") ? a*ones(N) : a*ones(N,2)

    # noise
    # L = length(cword)
    noise = Vector{Float64}(undef,N-twoZc)

    # received signal
    signal = zeros(N-twoZc)

    # bit-error
    biterror = Vector{Bool}(undef,N) 

    # Lq -> matrix of Nv messages (N x M)
    # Lr -> matrix of Nc messages (M x N)
    # if mode == "MKAY" the are matrices for bit = 0 and bit = 1
    Lq = (bptype != "MKAY") ? zeros(M,N) : zeros(M,N,2)

    Lr = (bptype != "MKAY") ? zeros(M,N) : zeros(M,N,2)
    
    # Set variables Lrn and signs depending on the BP type (used for dispatch)
    if bptype == "FAST" || bptype == "MKAY"
        Lrn = zeros(N)
        signs = nothing
    elseif bptype == "ALTN" || bptype == "TABL" 
        Lrn = zeros(N)
        signs = zeros(Bool,N)
    elseif bptype == "MSUM"
        Lrn = nothing
        signs = zeros(Bool,N)
    else
        Lrn = nothing
        signs = nothing
    end

    phi = (bptype == "TABL") ? PHI : nothing

    # Set other variables that depend on the mode
    newLr = RBP ? H*0.0 : nothing
    Factors = RBP ? 1.0*H  : nothing

    # RBP modes
    if mode == "RBP" || mode == "Genius-RBP"
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
        if mode == "Genius-RBP"
            global TOTALBITERROR = zeros(Int,num_edges)
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
            cword = NR_LDPC_encode(E_H,msg,Nr_ldpc_data)
        elseif protocol == "PEG"
            cword = gf2_mat_mult(G,msg)   
        elseif protocol == "WiMAX"
            cword = IEEE80216e_parity_bits(msg,Zf,E_H,E_M,E_N,E_K)
        end

        # 3) Modulation of the cword
        u = Float64.(2*cword .- 1)

        # 4) sum the noise to the modulated cword to produce the received signal
        received_signal!(signal,noise,stdev,u,rgn_noise,noisetest)

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
        
        # 9) init the Lq matrix
        init_Lq!(Lq,Lf,Nv)

        # 10) precalculate the residues for RBP
        if mode == "RBP" || mode == "List-RBP" || mode == "Genius-RBP"
            calc_all_residues!(Lq,Lr,Nc,Lrn,signs,phi,newLr,Factors,rbpmatrix,
                residues,coords,listsizes,relative)
        elseif mode == "NW-RBP"
            for ci in 1:M
                # calculate the new check to node messages
                Nci = Nc[ci]
                maxresidue = 0.0
                for vj in Nci
                    li = LinearIndices(Lr)[ci,vj]
                    newlr = calc_Lr(Nci,ci,vj,Lq)
                    newLr[li] = newlr
                    Lr[li] = newlr
                    residue = calc_residue(newlr,0.0,Factors[ci])
                    if residue > maxresidue
                        maxresidue = residue
                    end
                end
                residues[ci] = maxresidue
            end
        elseif mode == "VN-RBP"
            for m in 1:M 
                # calculate the new check to node messages
                update_Lr!(newLr,Lq,m,Nc[m],Lrn,signs,phi)
            end
            for n in 1:N
                residue = 0.0
                for m in Nv[n]
                    li = LinearIndices(newLr)[m,n]
                    residue += newLr[li]
                end
                if relative     
                    residue /= Lf[n]
                end
                residues[n] = abs(residue)
            end
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
                    Lrn,
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
                    Lrn,
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
                    Lrn,
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
            elseif mode == "Genius-RBP"
                bp_not_converged = genius_RBP!(
                    bitvector,
                    complete_cword,
                    Lq,
                    Lr,
                    Lf,
                    Nc,
                    Nv,
                    Lrn,
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
            elseif mode == "VN-RBP"
                bp_not_converged = VN_RBP!(
                    bitvector,
                    Lq,
                    Lr,
                    Lf,
                    Nc,
                    Nv,
                    Lrn,
                    decayfactor,
                    num_edges,
                    newLr,
                    Factors,
                    residues,
                    bp_not_converged
                    )
                # reset factors
                Factors .= 1.0

            elseif mode == "NW-RBP"
                bp_not_converged = NW_RBP!(
                    bitvector,
                    Lq,
                    Lr,
                    Lf,
                    Nc,
                    Nv,
                    Lrn,
                    decayfactor,
                    M,
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