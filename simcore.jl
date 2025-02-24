################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Core routine to estimate the LPCD performance (FER BER x SNR)

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
include("./RBP functions/find_local_maxresidue.jl")

function
    simcore(
        A::Integer,
        snr::Real,
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
    DECODED = zeros(Int,maxiter)
    decoded = Vector{Bool}(undef,maxiter)

    # bit error rate
    BER = zeros(Int,maxiter)
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
    all_max_res = (test && RBP) ? zeros(MAXRBP*num_edges) : nothing
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
    max_coords = (mode == "Local-RBP") ? zeros(Int,4) : nothing
    max_residue = (mode == "Local-RBP") ? ones(2)*0 : nothing
    max_coords_alt = (mode == "Local-RBP") ? zeros(Int,2) : nothing

    # mode = List-RBP
    if mode == "List-RBP" || mode == "Mod-List-RBP" || mode == "Random-List-RBP"
        listres = zeros(listsizes[1]+1)
        listm = zeros(Int,listsizes[1]+1)
        listn = zeros(Int,listsizes[1]+1)
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
        listres = nothing
        listm = nothing
        listn = nothing
        listres2 = nothing
        listm2 = nothing
        listn2 = nothing
        inlist = nothing
        listsizes[2] = 0
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
        bitvector .= true
        syndrome .= true
        decoded .= false
        for m in 1:M
            for n in cn2vn[m]
                Lr[m,n] = 0.0
            end
        end
        if mode == "List-RBP" || mode == "Mod-List-RBP"
            listres .= 0.0
            listm .= 0
            listn .= 0
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
            for m in eachindex(cn2vn)
                calc_residues!(addressinv,residues,Factors,Ms,nothing,Lq,Lrn,
                signs,phi,0,m,cn2vn)
            end
        elseif mode == "Local-RBP"
            for m in eachindex(cn2vn)
                find_local_maxresidue!(max_residue,Factors,Ms,nothing,Lq,Lrn,
                signs,phi,0,m,cn2vn,max_coords)
            end
            max_coords_alt[1] = max_coords[3]
            max_coords_alt[2] = max_coords[4]
        elseif mode == "List-RBP" || mode == "Mod-List-RBP" || mode == "Random-List-RBP"
            for m in eachindex(cn2vn)
                _ = calc_residues!(Factors,Ms,nothing,Lq,Lrn,
                signs,phi,0,m,cn2vn,listres,listm,listn,listres2,listm2,listn2,
                listsizes,0,inlist)
            end
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
                    max_coords,
                    max_coords_alt,
                )
                # reset factors
                resetfactors!(Factors,vn2cn)
            elseif mode == "List-RBP"
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
                    listres,
                    listm,
                    listn,
                    listres2,
                    listm2,
                    listn2,
                    inlist
                )
                # reset factors
                resetfactors!(Factors,vn2cn)            
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
                    listres,
                    listm,
                    listn,
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
                    listres,
                    listm,
                    listn,
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
                    println("Max LLR estimate errors: ")
                for j in eachindex(bitvector)
                    print(Int(bitvector[j] != cword[j]))
                    if j%80 == 0
                        println()
                    end
                end     
                println() 
                println("Syndrome: ")
                for j in eachindex(syndrome)
                    print(Int(syndrome[j]))
                    if j%80 == 0
                        println()
                    end
                end     
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
        @. BER += ber
        @. DECODED += decoded

    end

    if test
        return Lr, Lq, all_max_res
    else
        return DECODED, BER
    end

end

