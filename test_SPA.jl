################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Functions to test the SPA algorithms

function test_SPA(indices_m::Vector{Vector{Int64}}, 
                  indices_n::Vector{Vector{Int64}},
                  h::Matrix{Int64},
                  M::Int64,
                  N::Int64,
                  phi::Vector{Float64},
                  mode::String)

    σ = 0.8

    t = [1.3129, 2.6584, 0.7413, 2.1745, 0.5981, −0.8323, −0.3962, −1.7586,
    1.4905, 0.4084, −0.9290, 1.0765]

    f = zeros(N,2)
    ΔLf = zeros(N)
    k = 1/(sqrt(2π)*σ)
    s = 2σ^2
    for i in eachindex(t)
        f[i,1] = k*exp(-(t[i]+1)^2/s)
        f[i,2] = k*exp(-(t[i]-1)^2/s)
        if mode == "TABLE"
            ΔLf[i] = -SIZE_per_RANGE*4t[i]/s
        else
            ΔLf[i] = -4t[i]/s
        end
    end

    normalize(f,N)

    Q = zeros(N,M,2)
    Lq = zeros(N,M)

    Q = init_q(Q,N,f,indices_m)
    Lq = llr_init_q(Lq, N, ΔLf, indices_m)


    S = -1
    i = 0

    R = zeros(M,N,2)
    Lr = zeros(M,N)
    d_llr = zeros(Int64, N)

    Lrn = zeros(N)

    sn = ones(Int,N)

    while S != 0 && i < MAX

        i += 1

        println("Iteration #", i)

        println("Conventional SPA:")

        δQ = Q[:,:,1]-Q[:,:,2]

        d = zeros(Int, N)

        R = simple_horizontal_update(R, M, δQ, indices_n)
        Q, d = vertical_update_and_MAP(d, N,R,Q,f, indices_m)

        syndrome = h*d .% 2

        println("MAP estimate syndrome: ", syndrome)

        println("LLR SPA:")

        if mode == "TNH"
            # tanh SPA
            Lr = llr_horizontal_update(M,Lr,Lq,indices_n)
        elseif mode == "APP"
            # approximate SPA
            Lr = llr_horizontal_update(M,Lr,Lq,indices_n,sn)
        elseif mode == "ALT"
            # alternative SPA
            Lr = llr_horizontal_update(M,Lr,Lq,indices_n,sn,Lrn)
        elseif mode == "TAB"
            # lookup-table SPA
            Lr = llr_horizontal_update(M,Lr,Lq,indices_n,sn,Lrn,phi)
        else
            println("ERROR")
        end

        Lq, d_llr = llr_vertical_update_and_MAP(Lq,Lr,d_llr,N,ΔLf,indices_m)

        syndrome_llr = h*d_llr .% 2   

        println("LLR MAP estimate syndrome: ", syndrome_llr)
        
        S = sum(syndrome)

    end

    return R, Lr, Q, Lq
end