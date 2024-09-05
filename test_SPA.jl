################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Functions to test the SPA algorithms

function test_SPA(indices_m::Vector{Vector{Int64}}, 
                  indices_n::Vector{Vector{Int64}},
                  phi::Vector{Float64},
                  mode::String)

    M = length(indices_n)
    N = length(indices_m)

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

    init_q!(Q,N,f,indices_m)
    llr_init_q!(Lq,ΔLf,indices_m)


    S = -1
    i = 0

    R = zeros(M,N,2)
    Lr = zeros(M,N)
    d = zeros(Int, N)
    d_llr = zeros(Int64, N)

    Lrn = zeros(N)

    sn = ones(Int,N)

    syndrome = zeros(Int,M)
    syndrome_llr = zeros(Int,M)

    while S != 0 && i < MAX

        i += 1

        println("Iteration #", i)

        ### Conventional simplified SPA

        δQ = Q[:,:,1]-Q[:,:,2]    

        simple_horizontal_update!(R,δQ,indices_n)
        vertical_update_and_MAP!(Q,d,R,f,indices_m)

        calc_syndrome!(syndrome,indices_n,d)
        
        ### LLR SPA

        if mode == "TNH"
            # tanh SPA
            llr_horizontal_update!(Lr,Lq,indices_n,nothing,nothing,nothing)
        elseif mode == "APP"
            # approximate SPA
            llr_horizontal_update!(Lr,Lq,indices_n,sn,nothing,nothing)
        elseif mode == "ALT"
            # alternative SPA
            llr_horizontal_update!(Lr,Lq,indices_n,sn,Lrn,nothing)
        elseif mode == "TAB"
            # lookup-table SPA
            llr_horizontal_update!(Lr,Lq,indices_n,sn,Lrn,phi)
        else
            throw(ArgumentError("invalid mode"))
        end

        llr_vertical_update_and_MAP!(Lq,d_llr,Lr,ΔLf,indices_m)

        calc_syndrome!(syndrome_llr,indices_n,d)

        println("MAP SPL estimate: ", d)

        println("MAP LLR estimate: ", d_llr)

        println("MAP SPL estimate syndrome: ", syndrome)

        println("LLR MAP estimate syndrome: ", syndrome_llr)
        
        S = sum(syndrome)

    end

    return R, Lr, Q, Lq
end