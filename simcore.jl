################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Core routine to estimate the LPCD performance (FER BER x SNR)

include("auxiliary_functions.jl")
include("lookupTable.jl")
include("BP.jl")
include("calc_Lf.jl")
include("./RBP functions/calc_residues.jl")
include("./RBP functions/find_local_maxresidue.jl")

function
    performance_simcore(
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
        decay::Union{AbstractFloat,Nothing},
        listsize1::Integer,
        listsize2::Integer,
        rgn_seed_noise::Integer,
        rng_seed_sample::Integer,
        rgn_seed_msg::Integer;
        test=false,
        printtest=false
    )

    if mode == "RBP" || mode == "Local-RBP" || mode == "List-RBP" 
        supermode = "RBP"
    elseif mode == "LBP" || mode == "iLBP"
        supermode = "LBP"  
    else
        supermode = "Flooding"
    end
    
################################## CONSTANTS ###################################

    Zc = nr_ldpc_data.Zc

    # transform snr in standard deviations
    variance = 1 ./ (exp10.(snr/10))
    stdev = sqrt.(variance)
    
    # list of checks and variables nodes
    vn2cn  = make_vn2cn_list(H)
    cn2vn  = make_cn2vn_list(H)

    ############################# PREALLOCATIONS ###############################

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

    Ms = (supermode == "RBP") ? H*0.0 : nothing
    
    # Set variables Lrn and signs depending on the BP type (also used for dispatch)
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

   # Set other variables that depend on the mode

    visited_vns = (mode == "iLBP") ? zeros(Bool,N) : nothing

    if mode == "iLBP" || supermode == "RBP"       
        Ldn = zeros(N)
    else
        Ldn = nothing
    end

    phi = (bptype == "TABL") ? lookupTable() : nothing

    Factors, num_edges = (supermode == "RBP") ? (1.0*H, sum(H)) : (nothing,nothing)

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

    if mode == "List-RBP"
        listres1 = zeros(listsize1+1)
        listm1 = zeros(Int,listsize1+1)
        listn1 = zeros(Int,listsize1+1)
        inlist = Matrix(false*H)
        if listsize2 != 0
            listres2 = zeros(listsize2+1)
            listm2 = zeros(Int,listsize2+1)
            listn2 = zeros(Int,listsize2+1)
        else
            listres2 = nothing
            listm2 = nothing
            listn2 = nothing
        end
    else        
        if mode == "RBP" || mode == "Local-RBP"
            listsize1 = 1
            listres1 = zeros(1)
            listm1 = zeros(Int,1)
            listn1 = zeros(Int,1)
        else
            listres1 = nothing
            listm1 = nothing
            listn1 = nothing
        end
        listres2 = nothing
        listm2 = nothing
        listn2 = nothing
        inlist = nothing
        listsize2 = 0
    end 

################################## MAIN LOOP ###################################
    
    rng_noise = Xoshiro(rgn_seed_noise)

    rgn_msg = Xoshiro(rgn_seed_msg)

    rng_sample = (supermode == "RBP") ? Xoshiro(rng_seed_sample) : nothing  

    if mode == "Local-RBP"
        supermode = "Local-RBP"
    end
    maxcoords = zeros(Int,4)
    maxresidues = ones(2)*0
    maxcoords_alt = zeros(Int,2)

    @fastmath @inbounds for j in 1:trials

        rand!(rgn_msg,msg,Bool)
        if LDPC == 1
            cword = NR_LDPC_encode(E_H, msg, nr_ldpc_data)
        elseif LDPC == 3
            cword = IEEE80216e_parity_bits(msg,Zf,E_H)
        end

        if test && printtest
            println("Realization #$j:")
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
        u = Float64.(2*cword .- 1)
        if Zc > 0
            cword = [msg[1:2*Zc]; cword]
        end

        if mode == "List-RBP"
            listres1 .*= 0.0
            listm1 .*= 0
            listn1 .*= 0
            if listsize2 != 0
                listres2 .*= 0.0
                listm2 .*= 0
                listn2 .*= 0
            end
            for m in 1:M
                for n in cn2vn[m]
                    inlist[m,n] = false
                end
            end
        end

        received_signal!(signal,noise,stdev,u,rng_noise)

        bitvector .= true

        syndrome .= true

        decoded .= false
        
        # reinitialize matrices Lr and Lq
        for m in 1:M
            for n in cn2vn[m]
                Lr[m,n] *= 0.0
                Lq[n,m] *= 0.0
            end
        end

        # init the llr priors
        calc_Lf!(view(Lf,2*Zc+1:N),signal,variance)
        if bptype == "TABL"
            # scale for table
            Lf .*= SIZE_per_RANGE
        end
        init_Lq!(Lq,Lf,vn2cn)

        if supermode == "RBP"
            for m in eachindex(cn2vn)
                calc_residues!(addressinv,residues,nothing,Ms,nothing,Lq,Lrn,signs,phi,0,m,
                cn2vn,listres1,listm1,listn1,nothing,nothing,nothing,listsize1,0,
                inlist)
            end
        elseif supermode == "Local-RBP"
            for m in eachindex(cn2vn)
                find_local_maxresidue!(maxresidues,nothing,Ms,nothing,Lq,Lrn,signs,phi,0,m,
                cn2vn,maxcoords)
            end
            maxcoords_alt[1] = maxcoords[3]
            maxcoords_alt[2] = maxcoords[4]
        end
            
        # SPA routine
        BP!(address,
            addressinv,
            supermode,
            stop,
            test,
            maxiter,
            syndrome,
            bitvector,
            cword,
            biterror,
            ber,
            decoded,
            Lf,
            Lq,
            Lr,
            Ms,
            cn2vn,
            vn2cn,
            Lrn,
            signs,
            phi,
            printtest,
            residues,
            Factors,
            decay,
            num_edges,
            Ldn,
            visited_vns,
            rng_sample,
            listsize1,
            listsize2,
            listres1,
            listm1,
            listn1,
            listres2,
            listm2,
            listn2,
            inlist,
            maxresidues,
            maxcoords,
            maxcoords_alt
            )                


        # bit error rate
        @. BER += ber
        @. DECODED += decoded

    end

    if test
        return Lr, Lq
    else
        return DECODED, BER
    end

end

