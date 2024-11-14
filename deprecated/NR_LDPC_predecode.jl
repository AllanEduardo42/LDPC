################################################################################
# Allan Eduardo Feitosa
# 13 NOV 2024
# Function to prepare a message to be decoded by the NR-LDPC

include("rate_matching.jl")

function
    NR_LDPC_predecode(
        g::Vector{Bool},
        nr_ldpc::NR_LDPC;
        Q_m = 1
    )

    E_r = nr_ldpc.E_r
    C = nr_ldpc.C
    N = nr_ldpc.N
    N_cb = nr_ldpc.N_cb
    k0 = nr_ldpc.k0
    K = nr_ldpc.K
    K_prime = nr_ldpc.K_prime
    Zc = nr_ldpc.Zc
    
    # Code concatenation

    f = inv_code_concatenation(C,E_r,g)

    # bit interleaving

    e = inv_bit_interleaving(C,E_r,Q_m,f)

    # bit selection
    range = (K_prime-2*Zc+1):(K-2*Zc)
    d = inv_rate_matching(E_r,N,C,e,N_cb,k0,range)

    return d

end