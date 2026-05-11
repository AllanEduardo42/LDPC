function print_simulation_details(
    errors::Int,
    maxiter::Int,
    ebn0::Union{Float64,Vector{Float64}}
)
    str = """
********************************************************************************

                           LDPC PERFORMANCE SIMULATION                          

********************************************************************************

############################### LDPC parameters ################################

LDPC Protocol: """

    if PROTOCOL == "5GNR"
        str *= "NR-LDPC (5G), Zc = $(NR_LDPC_DATA.Zc), BG = $(NR_LDPC_DATA.bg), set index = $(NR_LDPC_DATA.iLS)"
    elseif PROTOCOL == "PEG"
        str *= "PEG"
    elseif PROTOCOL == "WiMAX"
        str *= "IEEE80216e (WiMAX)"
    end
    str *= "\nParity Check Matrix: $MM x $NN\n"

    println(str)
    if SAVE
        println(FILE,str)
    end

    display(sparse(HH))

    str = """

Graph girth = $GIRTH
Rate = $(RATE[1])/$(RATE[2])
Code length = $CODE_LENGTH    

############################ Simulation parameters #############################

Maximum number of frame errors: $errors
Maximum number of iterations: $maxiter
Number of threads (multithreading): $NTHREADS
Eb/No (dB): $ebn0
"""
    if RAYL
        str*= """Rayleigh fading channel model activated

        """
    else
        str*= """
        
        """
    end
    if TEST
        str *= """######################### Starting simulation (Testing) ########################
        
        """
    else
        str *= """############################# Starting simulation ##############################

        """
    end
    
    print(str)
    if SAVE
        print(FILE,str)
    end

end