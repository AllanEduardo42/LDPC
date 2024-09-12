################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Function to test the SPA algorithms

function 
    test_SPA(
        c::Vector{Bool},
        indices_col::Vector{Vector{Int64}}, 
        indices_row::Vector{Vector{Int64}},
        t::Vector{Float64},
        σ::Float64,
        phi::Vector{Float64},
        mode::String,
        printing::Bool
    )

    M = length(indices_row)
    N = length(indices_col)

    f = zeros(N,2)
    ΔLf = zeros(N)
    k = 1/(sqrt(2π)*σ)
    s = σ^2
    for i in eachindex(t)
        f[i,1] = k*exp(-(t[i]+1)^2/(2*s))
        f[i,2] = k*exp(-(t[i]-1)^2/(2*s))
        if mode == "TAB"
            ΔLf[i] = -2*SIZE_per_RANGE*t[i]/s
        else
            ΔLf[i] = -2t[i]/s
        end
    end

    normalize!(f)

    q = zeros(N,M,2)
    Lq = zeros(N,M)

    init_q!(q,f,indices_col)
    llr_init_q!(Lq,ΔLf,indices_col)

    r = zeros(M,N,2)
    Lr = zeros(M,N)
    d = zeros(Bool, N)
    d_llr = zeros(Bool, N)

    Lrn = zeros(N)

    sn = ones(Bool,N)

    syndrome = ones(Bool,M)
    syndrome_llr = zeros(Bool,M)

    index = MAX
    FIRST = true
    DECODED = false
    i = 1
    converged = false
    while i <= MAX && !converged

        ### Conventional simplified SPA

        δQ = q[:,:,1]-q[:,:,2]    

        simple_horizontal_update!(
            r,
            δQ,
            indices_row
        )
        vertical_update_and_MAP!(
            q,
            d,
            r,
            f,
            indices_col
        )

        calc_syndrome!(
            syndrome,
            d,
            indices_row
        )
        
        ### LLR SPA

        if mode == "TNH"
            # tanh SPA
            llr_horizontal_update!(
                Lr,
                Lq,
                indices_row,
                nothing,
                nothing,
                nothing
            )
        elseif mode == "ALT"
            # alternative SPA
            llr_horizontal_update!(
                Lr,
                Lq,
                indices_row,
                sn,
                Lrn,
                nothing
            )
        elseif mode == "TAB"
            # lookup-table SPA
            llr_horizontal_update!(
                Lr,
                Lq,
                indices_row,
                sn,
                Lrn,
                phi
            )
        elseif mode == "MIN"
            # approximate SPA
            llr_horizontal_update!(
                Lr,
                Lq,
                indices_row,
                sn,
                nothing,
                nothing
            )
        else
            throw(
                ArgumentError(
                    "$mode is not a valid mode"
                )
            )
        end

        converged = llr_vertical_update_and_MAP!(
            Lq,
            d_llr,
            Lr,
            ΔLf,
            indices_col
        )
        println(converged)

        calc_syndrome!(
            syndrome_llr,
            d_llr,
            indices_row
        )

        if printing

            println("Iteration #$i")

            println("MAP SIMPLE estimate: $d")

            println("MAP Δ-LLRs estimate: $d_llr")

            println("MAP SIMPLE syndrome: $syndrome")

            println("LLR Δ-LLRs syndrome: $syndrome_llr")

        end

        if FIRST && iszero(syndrome)
            FIRST = false
            index = i
            if d == c
                DECODED = true
            end
        end

        i += 1

    end

    return r, Lr, q, Lq, index, DECODED
end
