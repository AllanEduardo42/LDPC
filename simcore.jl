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
include("Local_RBP.jl")
include("List-RBP.jl")
include("Mod-List-RBP.jl")
include("Random-List-RBP.jl")
include("./RBP functions/calc_residues.jl")

function
    simcore(
        A::Integer,
        snr::AbstractFloat,
        H::BitMatrix,
        E_H::Matrix{<:Integer},
        LDPC::Integer,
        Zf::Integer,
        nr_ldpc_data::NR_LDPC_DATA,
        mode::String,
        bptype::String,
        trials::Integer,
        maxiter::Integer,
        stop::Bool,
        decayfactor::AbstractFloat,
        listsizes::Vector{<:Integer},
        rgn_seed_noise::Integer,
        rng_seed_sample::Integer,
        rgn_seed_msg::Integer;
        test=false,
        printtest=false
    )

    if mode == "RBP" || mode == "Local-RBP" || mode == "List-RBP" || 
       mode == "Mod-List-RBP" || mode == "Random-List-RBP"
        RBP = true
    else
        RBP = false
    end
    
################################## CONSTANTS ###################################

    M,N = size(H)

    # constant Zc of NR/5G
    Zc = nr_ldpc_data.Zc

    # number of edges in the graph
    num_edges = sum(H)

    # transform snr in standard deviations
    variance = 1 ./ (exp10.(snr/10))
    stdev = sqrt.(variance)
    
    # list of checks and variables nodes
    vn2cn  = make_vn2cn_list(H)
    cn2vn  = make_cn2vn_list(H)

    # Set the random seeds
    rng_noise = Xoshiro(rgn_seed_noise)
    rgn_msg = Xoshiro(rgn_seed_msg)
    rng_rbp = RBP ? Xoshiro(rng_seed_sample) : nothing  

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
    Lf = (bptype != "MKAY") ? zeros(N) : zeros(N,2)

    # noise
    # L = length(cword)
    noise = Vector{Float64}(undef,N-2*Zc)

    # received signal
    signal = zeros(N-2*Zc)

    # bit-error
    biterror = Vector{Bool}(undef,N) 

    # Lq -> matrix of vn2cn messages (N x M)
    # Lr -> matrix of cn2vn messages (M x N)
    # if mode == "MKAY" the are matrices for bit = 0 and bit = 1
    Lq = (bptype != "MKAY") ? zeros(N,M) : zeros(N,M,2)

    Lr = (bptype != "MKAY") ? zeros(M,N) : zeros(M,N,2)
    
    # Set variables Lrn and signs depending on the BP type (used for dispatch)
    if bptype == "FAST"
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

    phi = (bptype == "TABL") ? lookupTable() : nothing

    # Set other variables that depend on the mode

    # supermode = RBP
    Ms = RBP ? H*0.0 : nothing
    Factors  = RBP ? 1.0*H  : nothing
    all_max_res = (test && RBP) ? zeros(maxiter*num_edges) : nothing
    all_max_res_alt = (test && RBP) ? zeros(num_edges) : nothing

    # mode = RBP
    if mode == "RBP"
        residues = zeros(num_edges)
        address = zeros(Int,2,num_edges)
        addressinv = 0*H
        e = 0
        for m in axes(H,1)
            for n in cn2vn[m]
                e += 1
                address[1,e] = m
                address[2,e] = n
                addressinv[m,n] = e
            end
        end
    else
        residues = nothing
        address = nothing
        addressinv = nothing
    end

    # mode = Local-RBP
    max_coords = (mode == "Local-RBP") ? zeros(Int,6) : nothing
    max_residue = (mode == "Local-RBP") ? zeros(3) : nothing

    # mode = List-RBP
    if mode == "List-RBP" || mode == "Mod-List-RBP" || mode == "Random-List-RBP"
        listres1 = zeros(listsizes[1]+1)
        listm1 = zeros(Int,listsizes[1]+1)
        listn1 = zeros(Int,listsizes[1]+1)
        inlist = Matrix(false*H)
        if listsizes[2] > 0
            if listsizes[2] == 1
                listres2 = zeros(listsizes[1]+1)
                listm2 = zeros(Int,listsizes[1]+1)
                listn2 = zeros(Int,listsizes[1]+1)
            else
                listres2 = zeros(listsizes[2]+1)
                listm2 = zeros(Int,listsizes[2]+1)
                listn2 = zeros(Int,listsizes[2]+1)
            end
        else
            listres2 = nothing
            listm2 = nothing
            listn2 = nothing
        end
    else        
        listres1 = nothing
        listm1 = nothing
        listn1 = nothing
        listres2 = nothing
        listm2 = nothing
        listn2 = nothing
        inlist = nothing
    end

################################## MAIN LOOP ###################################

    @inbounds for j in 1:trials

        # generate the random message
        rand!(rgn_msg,msg,Bool)

        # generate the cword
        if LDPC == 1
            cword = NR_LDPC_encode(E_H, msg, nr_ldpc_data)
        elseif LDPC == 3
            cword = IEEE80216e_parity_bits(msg,Zf,E_H)
        end

        # Modulation of the cword
        u = Float64.(2*cword .- 1)

        # Include the punctured bits in the cword for biterror calculation
        if Zc > 0
            cword = [msg[1:2*Zc]; cword]
        end

        # sum the noise to the modulated cword to produce the received signal
        received_signal!(signal,noise,stdev,u,rng_noise)

        # print info in test mode
        if test && printtest
            println("Trial #$j:")
            println()
            println("msg (L = $(length(msg))):")
            for i in eachindex(msg)
                print(Int(msg[i]))
                if i%80 == 0
                    println()
                end
            end
            println()
            println()
            println("cword (L = $(length(cword))):")
            for i in eachindex(cword)
                print(Int(cword[i]))
                if i%80 == 0
                    println()
                end
            end
            println()
            println()
        end       

        # reset the simulation variables
        bitvector .= false
        syndrome .= true
        decoded .= false
        for m in 1:M
            for n in cn2vn[m]
                Lr[m,n] = 0.0
            end
        end
        if mode == "List-RBP" || mode == "Mod-List-RBP"
            listres1 .= 0.0
            listm1 .= 0
            listn1 .= 0
            if listsizes[2] != 0
                listres2 .= 0.0
                listm2 .= 0
                listn2 .= 0
            end
            for m in 1:M
                for n in cn2vn[m]
                    inlist[m,n] = false
                end
            end
        end

        # init the llr priors
        calc_Lf!(view(Lf,2*Zc+1:N),signal,variance)
        if bptype == "TABL"
            # scale for table
            Lf .*= SIZE_per_RANGE
        end
        # init the Lq matrix
        init_Lq!(Lq,Lf,vn2cn)

        # precalculate the residues in RBP
        if mode == "RBP"
            rbp_not_converged = true
            calc_residues!(addressinv,residues,Factors,Ms,Lr,Lq,Lrn,signs,phi,0,
                0,cn2vn,eachindex(cn2vn))
        elseif mode == "Local-RBP"
            rbp_not_converged = true
            calc_residues!(max_residue,max_coords,Factors,Ms,Lr,Lq,Lrn,signs,
                phi,0,0,cn2vn,eachindex(cn2vn))
        elseif mode == "List-RBP"
            rbp_not_converged = true
            calc_residues!(Factors,Ms,Lr,Lq,Lrn,signs,phi,0,0,cn2vn,
                eachindex(cn2vn),listres1,listm1,listn1,listres2,listm2,listn2,
                [listsizes[1], 0],inlist)
        end
        
        # BP routine
        
        for iter in 1:maxiter

            if test && printtest  
                println("### Iteration #$iter ###")
            end
    
            if mode == "Flooding"
                flooding!(
                    bitvector,
                    Lq,
                    Lr,
                    Lf,
                    cn2vn,
                    vn2cn,
                    Lrn,
                    signs,
                    phi)  
            elseif mode == "LBP"
                LBP!(
                    bitvector,
                    Lq,
                    Lr,
                    Lf,
                    cn2vn,
                    vn2cn,
                    Lrn)
            elseif mode == "RBP"
                RBP!(
                    bitvector,
                    Lq,
                    Lr,
                    Lf,
                    cn2vn,
                    vn2cn,
                    Lrn,
                    signs,
                    phi,
                    decayfactor,
                    num_edges,
                    Ms,
                    Factors,
                    all_max_res_alt,
                    test,
                    address,
                    addressinv,
                    residues,
                    )
                # reset factors
                resetfactors!(Factors,vn2cn)
            elseif mode == "Local-RBP"
                if rbp_not_converged
                    local_RBP!(
                        bitvector,
                        Lq,
                        Lr,
                        Lf,
                        cn2vn,
                        vn2cn,
                        Lrn,
                        signs,
                        phi,
                        decayfactor,
                        num_edges,
                        Ms,
                        Factors,
                        all_max_res_alt,
                        test,
                        rng_rbp,
                        max_residue,
                        max_coords
                    )
                    # reset factors
                    resetfactors!(Factors,vn2cn)
                    if max_residue[1] == 0.0
                        calc_residues!(max_residue,max_coords,Factors,Ms,Lr,Lq,
                            Lrn,signs,phi,0,0,cn2vn,eachindex(cn2vn))
                        if max_residue[1] == 0.0
                            rbp_not_converged = false
                        end
                    end
                end              
            elseif mode == "List-RBP"
                if rbp_not_converged
                    list_RBP!(
                        bitvector,
                        Lq,
                        Lr,
                        Lf,
                        cn2vn,
                        vn2cn,
                        Lrn,
                        signs,
                        phi,
                        decayfactor,
                        num_edges,
                        Ms,
                        Factors,
                        all_max_res_alt,
                        test,
                        rng_rbp,
                        listsizes,
                        listres1,
                        listm1,
                        listn1,
                        listres2,
                        listm2,
                        listn2,
                        inlist
                    )
                    # reset factors
                    resetfactors!(Factors,vn2cn)
                    if listres1[1] == 0.0
                        calc_residues!(Factors,Ms,Lr,Lq,Lrn,signs,phi,0,0,cn2vn,
                            eachindex(cn2vn),listres1,listm1,listn1,listres2,listm2,
                            listn2,[listsizes[1], 0],inlist)
                        if listres1[1] == 0.0
                            rbp_not_converged = false
                        end
                    end
                end       
            elseif mode == "Mod-List-RBP"
                mod_list_RBP!(
                    bitvector,
                    Lq,
                    Lr,
                    Lf,
                    cn2vn,
                    vn2cn,
                    Lrn,
                    signs,
                    phi,
                    decayfactor,
                    num_edges,
                    Ms,
                    Factors,
                    all_max_res_alt,
                    test,
                    rng_rbp,
                    listsizes,
                    listres1,
                    listm1,
                    listn1,
                    listres2,
                    listm2,
                    listn2,
                    inlist,
                    syndrome
                )
                # reset factors
                resetfactors!(Factors,vn2cn)
            elseif mode == "Random-List-RBP"
                random_list_RBP!(
                    bitvector,
                    Lq,
                    Lr,
                    Lf,
                    cn2vn,
                    vn2cn,
                    Lrn,
                    signs,
                    phi,
                    decayfactor,
                    num_edges,
                    Ms,
                    Factors,
                    all_max_res_alt,
                    test,
                    rng_rbp,
                    listsizes,
                    listres1,
                    listm1,
                    listn1,
                    listres2,
                    listm2,
                    listn2,
                    inlist
                )
                # reset factors
                resetfactors!(Factors,vn2cn)
            end
    
            calc_syndrome!(syndrome,bitvector,cn2vn)
    
            if test && printtest
    
                if RBP
                    all_max_res[(1 + (iter-1)*num_edges):(iter*num_edges)] = all_max_res_alt
                end
                    biterror .= (bitvector .≠ cword)
                    println("Bit error:")
                for j in eachindex(bitvector)
                    print(Int(biterror[j]))
                    if j%80 == 0
                        println()
                    end
                end
                println() 
                println("Bit error rate: $(sum(biterror))/$N")     
                println() 
                println("Syndrome: ")
                for j in eachindex(syndrome)
                    print(Int(syndrome[j]))
                    if j%80 == 0
                        println()
                    end
                end     
                println()
                println("Syndrome rate: $(sum(syndrome))/$M")     
                println() 
                if iszero(syndrome) && stop
                    break
                end
                println()     
            else
                if iszero(syndrome)
                    if bitvector == cword
                        @inbounds decoded[iter] = true
                    end
                    if stop
                        if iter < maxiter
                            @inbounds decoded[iter+1:end] .= decoded[iter]
                        end
                        break
                    end
                end
                biterror .= (bitvector .≠ cword)
                @inbounds ber[iter] = sum(biterror)
            end
        end

        # bit error rate
        @. sum_ber += ber
        @. sum_decoded += decoded

    end

    if test
        return Lr, Lq, all_max_res
    else
        return sum_decoded, sum_ber
    end

end

