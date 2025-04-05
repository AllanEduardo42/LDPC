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
include("VLBP.jl")
include("RBP.jl")
include("Local_RBP.jl")
include("List-RBP.jl")
# include("Mod-List-RBP.jl")
# include("Random-List-RBP.jl")
include("NRBP.jl")
include("./RBP functions/calc_all_residues.jl")

function
    simcore(
        A::Integer,
        snr::AbstractFloat,
        H::Matrix{Bool},
        G::Union{Nothing,Matrix{Bool}},
        M::Integer,
        N::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        E_H::Union{Nothing,Matrix{<:Integer}},
        LDPC::Integer,
        Zf::Integer,
        Nr_ldpc_data::Union{Nothing,nr_ldpc_data},
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
        printtest=false,
        msgtest=nothing,
        noisetest=nothing        
    )

    if mode == "RBP" || mode == "Local-RBP" || mode == "List-RBP" || 
       mode == "Mod-List-RBP" || mode == "Random-List-RBP" || 
       mode == "NRBP"
        RBP = true
    else
        RBP = false
    end
    
################################## CONSTANTS ###################################

    # constant Zc of NR/5G
    if LDPC == 1
        Zc = Nr_ldpc_data.Zc
    else
        Zc = 0
    end

    if LDPC == 3
        E_M, E_N = size(E_H)
        E_K = E_N - E_M
    end

    # number of edges in the graph
    num_edges = sum(H)

    # transform snr in standard deviations
    variance = 1 ./ (exp10.(snr/10))
    stdev = sqrt.(variance)

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
    Lf = (bptype != "MKAY") ? 0.01*ones(N) : 0.5*ones(N,2)

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
    Ms = RBP ? H*0.0 : nothing
    Factors = RBP ? 1.0*H  : nothing

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
    elseif mode == "NRBP"
        residues = zeros(N)
        Factors = ones(N)
    end

    # mode = Local-RBP
    maxcoords = (mode == "Local-RBP") ? zeros(Int,6) : nothing
    max_residue = (mode == "Local-RBP") ? zeros(3) : nothing

    # mode = List-RBP
    if mode == "List-RBP" || mode == "Mod-List-RBP"||
        mode == "Random-List-RBP"

        listres1 = zeros(listsizes[1]+1)
        coords1 = zeros(Int,2,listsizes[1]+1)
        inlist = Matrix(false*H)
        if listsizes[2] == 1
            listres2 = zeros(listsizes[1]+1)
            coords2 = zeros(Int,2,listsizes[1]+1)
        else
            listres2 = zeros(listsizes[2]+1)
            coords2 = zeros(Int,2,listsizes[2]+1)
        end
    end

################################## MAIN LOOP ###################################

    @inbounds for j in 1:trials

        # 1) generate the random message
        generate_message!(msg,rgn_msg,msgtest)

        # 2) generate the cword
        if LDPC == 1
            cword = NR_LDPC_encode(E_H,msg,Nr_ldpc_data)
        elseif LDPC == 2
            cword = gf2_mat_mult(G,msg)   
        elseif LDPC == 3
            cword = IEEE80216e_parity_bits(msg,Zf,E_H,E_M,E_N,E_K)
        end

        # 3) Modulation of the cword
        u = Float64.(2*cword .- 1)

        # 4) Include the punctured bits in the cword for biterror calculation
        if Zc > 0
            cword = [msg[1:2*Zc]; cword]
        end

        # 5) sum the noise to the modulated cword to produce the received signal
        received_signal!(signal,noise,stdev,u,rng_noise,noisetest)

        # 6) print info in test mode
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

        # 7) reset the simulation variables
        bitvector .= false
        syndrome .= true
        decoded .= false
        resetmatrix!(Lr,vn2cn,0.0)
        if mode == "List-RBP" || mode == "Mod-List-RBP"
            listres1 .= 0.0
            coords1 .= 0
            listres2 .= 0.0
            coords2 .= 0
            resetmatrix!(inlist,vn2cn,false)
        end

        # 8) init the llr priors
        if bptype == "MKAY"
            calc_Lf!(view(Lf,2*Zc+1:N,:),signal,variance)
        else
            calc_Lf!(view(Lf,2*Zc+1:N),signal,variance)
            if bptype == "TABL"
                # scale for table
                Lf .*= SIZE_PER_RANGE
            end
        end
        
        # 9) init the Lq matrix
        init_Lq!(Lq,Lf,vn2cn)

        # 10) precalculate the residues in RBP
        if mode == "RBP"
            calc_all_residues!(Lq,Lr,cn2vn,Lrn,signs,phi,Ms,Factors,addressinv,
                residues,M)
        # elseif mode == "NRBP"
        #     for m in 1:M 
        #         # calculate the new check to node messages
        #         update_Lr!(Ms,Lq,m,cn2vn,Lrn,signs,phi)
        #     end
        #     for n in 1:N
        #         for m in vn2cn[n]
        #             li = LinearIndices(Ms)[m,n]        
        #             residues[n] += Ms[li]
        #         end
        #         residues[n] = abs(residues[n]/Lf[n])
        #     end
        end

        # display(sort(residues,rev=true))
        
        # BP routine
        rbp_not_converged = true

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
            elseif mode == "VLBP"
                VLBP!(
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
                    address,
                    addressinv,
                    residues
                    )
                # reset factors
                resetmatrix!(Factors,vn2cn,1.0)
            elseif mode == "NRBP"
                NRBP!(
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
                    residues
                    )
                # reset factors
                Factors .= 1.0
            elseif mode == "Local-RBP"
                if rbp_not_converged && max_residue[1] == 0.0
                    calc_all_residues_local!(Lq,Lr,cn2vn,Lrn,signs,phi,Ms,Factors,
                        max_residue,maxcoords,M)
                    max_residue[3] = max_residue[2]
                    maxcoords[6] = maxcoords[4]
                    maxcoords[5] = maxcoords[3]
                    if max_residue[1] == 0.0
                        rbp_not_converged = false
                    end
                end
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
                        rng_rbp,
                        max_residue,
                        maxcoords
                    )
                    # reset factors
                    resetmatrix!(Factors,vn2cn,1.0)
                end              
            elseif mode == "List-RBP"
                if rbp_not_converged && listres1[1] == 0.0
                    calc_all_residues_list!(Lq,Lr,cn2vn,Lrn,signs,phi,Ms,Factors,
                        listsizes,listres1,coords1,inlist,M)
                    if listres1[1] == 0.0
                        rbp_not_converged = false
                    end
                end
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
                        rng_rbp,
                        listsizes,
                        listres1,
                        coords1,
                        listres2,
                        coords2,
                        inlist
                    )
                    # reset factors
                    resetmatrix!(Factors,vn2cn,1.0)
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
                    rng_rbp,
                    listsizes,
                    listres1,
                    coords1,
                    listn1,
                    listres2,
                    coords2,
                    listn2,
                    inlist,
                    syndrome
                )
                # reset factors
                resetmatrix!(Factors,vn2cn,1.0)
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
                    rng_rbp,
                    listsizes,
                    listres1,
                    coords1,
                    listn1,
                    listres2,
                    coords2,
                    listn2,
                    inlist
                )
                # reset factors
                resetmatrix!(Factors,vn2cn,1.0)
            end
    
            calc_syndrome!(syndrome,bitvector,cn2vn)
    
            if test && printtest
    
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
            if iter > 30
                thres = 0.0
            end
        end

        # bit error rate
        @. sum_ber += ber
        @. sum_decoded += decoded

    end

    if test
        if bptype == "MKAY"
            retr = zeros(M,N)
            retq = zeros(M,N)
            for m in eachindex(cn2vn)
                for n in cn2vn[m]
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

function generate_message!(
    msg::Vector{Bool},
    rgn_msg::AbstractRNG,
    ::Nothing
)

    rand!(rgn_msg,msg,Bool)

end

function generate_message!(
    msg::Vector{Bool},
    ::AbstractRNG,
    msgtest::Vector{Bool}
)

    msg .= msgtest

end