################################################################################
# Allan Eduardo Feitosa
# 2 set 2025
# C2V Update

####################### SPA USING FAST HYPERBOLIC TANGENT ######################
function calc_C2V(
    Nci::Vector{Int},
    ci::Int,
    vj::Int,    
    V2C::Matrix{Float64},
    ::Nothing
)
    
    @inbounds begin
        prod_v2c = 1.0
        for vb in Nci
            if vb ≠ vj
                v2c = V2C[ci,vb]
                if v2c == 0.0
                    return 0.0
                else
                    prod_v2c *= v2c
                end
            end
        end
        if abs(prod_v2c) < 1.0
            return 2*atanh(prod_v2c)
        elseif signbit(prod_v2c)
            return MINC2V
        else
            return MAXC2V
        end
    end
end

function calc_C2V_no_opt(
    Nci::Vector{Int},
    ci::Int,
    vj::Int,    
    V2C::Matrix{Float64},
    ::Union{Nothing,Float64}
)

    @fastmath @inbounds begin
        prod_v2c = 1.0
        for vb in Nci
            if vb ≠ vj
                v2c = V2C[ci,vb]
                if v2c == 0.0
                    return 0.0
                else
                    prod_v2c *= tanh(0.5*v2c)
                end
            end
        end
        if abs(prod_v2c) < 1.0
            return 2*atanh(prod_v2c)
        elseif signbit(prod_v2c)
            return MINC2V
        else
            return MAXC2V
        end
    end
end

############################### SPA USING MIN-SUM ##############################
function calc_C2V(
    Nci::Vector{Int},
    ci::Int,
    vj::Int,    
    V2C::Matrix{Float64},
    msum_factor::Float64
)

    @fastmath @inbounds begin
        s = false
        minL = MAXC2V
        for vb in Nci
            if vb ≠ vj
                v2c = V2C[ci,vb]
                sig = signbit(v2c)
                s ⊻= sig
                β = abs(v2c)
                if β < minL
                    minL = β
                end
            end
        end
        if s
            return -minL*msum_factor
        else
            return minL*msum_factor
        end
    end
end





