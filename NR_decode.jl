using SparseArrays

include("GF2_poly.jl")
include("make_parity_check_matrix.jl")
include("GF2_functions.jl")
include("rate_matching.jl")
include("encode.jl")

const R_LBRM = 2//3

function
    NR_LDPC_decode(
        g::Vector{Bool},
        E_r,
        C,
        k0,
        N_cb,
        N,
        K,
        K_prime,
        Zc;
        # R::Rational;
        # I_LBRM = 0,
        # TBS_LBRM = Inf,
        # rv = 0,
        # CBGTI = [],
        # N_L = 1,
        Q_m = 1
    )
    
    # Code concatenation

    f = inv_code_concatenation(C,E_r,g)

    # bit interleaving

    e = inv_bit_interleaving(C,E_r,Q_m,f)

    # bit selection
    d = inv_rate_matching(E_r,N,C,e,N_cb,k0,K,K_prime,Zc)

    # # LDPC decoder
    # for r = 1:C
        
    #     cw = [zeros(2*Z_c_,1); d_tilde[:,r]]
    #     c = cw[1:K]
    #     # cw(isnan(cw)) = inf;
    #     c_hat{r+1} = double(step(obj.hLDPCDecoder, cw));
    #     c_hat{r+1}(isnan(c_tilde)) = NaN;
    # end

    return d


end