################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Function to test the SPA algorithms

function 
    test_SPA(
        c::Vector{Bool},
        nodes2checks::Vector{Vector{T}} where {T<:Integer}, 
        checks2nodes::Vector{Vector{T}} where {T<:Integer},
        t::Vector{<:AbstractFloat},
        σ::AbstractFloat,
        phi::Vector{<:AbstractFloat},
        mode::String,
        printing::Bool
    )

    M = length(checks2nodes)
    N = length(nodes2checks)

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

    q = zeros(M,N,2)
    Lq = zeros(M,N)

    init_q!(q,f,nodes2checks)
    llr_init_q!(Lq,ΔLf,nodes2checks)

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
    while i <= MAX

        ### Conventional simplified SPA

        δQ = q[:,:,1]-q[:,:,2]    

        simple_horizontal_update!(
            r,
            δQ,
            checks2nodes
        )
        vertical_update_and_MAP!(
            q,
            d,
            r,
            f,
            nodes2checks
        )

        calc_syndrome!(
            syndrome,
            d,
            checks2nodes
        )
        
        ### LLR SPA

        if mode == "TNH"
            # tanh SPA
            llr_horizontal_update!(
                Lr,
                Lq,
                checks2nodes,
                nothing,
                nothing,
                nothing
            )
        elseif mode == "ALT"
            # alternative SPA
            llr_horizontal_update!(
                Lr,
                Lq,
                checks2nodes,
                sn,
                Lrn,
                nothing
            )
        elseif mode == "TAB"
            # lookup-table SPA
            llr_horizontal_update!(
                Lr,
                Lq,
                checks2nodes,
                sn,
                Lrn,
                phi
            )
        elseif mode == "MIN"
            # approximate SPA
            llr_horizontal_update!(
                Lr,
                Lq,
                checks2nodes,
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

        llr_vertical_update_and_MAP_crude!(
            Lq,
            d_llr,
            Lr,
            ΔLf,
            nodes2checks
        )

        calc_syndrome!(
            syndrome_llr,
            d_llr,
            checks2nodes
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
