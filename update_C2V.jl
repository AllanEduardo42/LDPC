################################################################################
# Allan Eduardo Feitosa
# 2 set 2025
# C2V Update
# There are 4 different methods for SPA: Mckay, tanh, minsum and table

####################### SPA USING FAST HYPERBOLIC TANGENT ######################
function 
    calc_C2V(
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

function 
    calc_C2V_no_opt(
        Nci::Vector{Int},
        ci::Int,
        vj::Int,    
        V2C::Matrix{Float64}
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
function 
    calc_C2V(
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

################ ALTERNATIVE TO HYPERBOLIC TANGENT USING TABLE ################
function 
    calc_ABCD!(
        V2C::Matrix{Float64},
        ci::Int,
        Nci::Vector{Int},
        signs::Vector{Bool},
        phi::Vector{Float64}
    )

    sum_c2v = 0.0
    s = false
    @fastmath @inbounds for vj in Nci
        v2c = V2C[ci,vj]
        sig = signbit(v2c)
        s ⊻= sig
        signs[vj] = sig
        sum_c2v += ϕ(abs(v2c),phi) 
    end

    return sum_c2v, s, nothing, nothing
end

function 
    calc_C2V(
        sum_c2v::Float64,   #A
        s::Bool,            #B          
        ::Nothing,          #C
        ::Nothing,          #D
        vj::Int,
        v2c::Float64,
        signs::Vector{Bool},
        phi::Vector{Float64}
    )

    @fastmath @inbounds begin
        x = abs(sum_c2v - v2c)
        y = signs[vj] ⊻ s
        return (1 - 2*y)*ϕ(x,phi)
    end
end

########################### SPA USING MKAY's METHOD ############################
function
    calc_δr(
        Nci::Vector{Int},
        vj::Int,
        δq::Vector{Float64}
    )
    
    @fastmath @inbounds begin
        δr = 1.0
        for vb in Nci
            if vb ≠ vj
                δr *= δq[vb]
            end
        end
        return δr
    end 

end

# pre-historic method
function
    calc_r(
        q::Array{Float64,3},
        ci::Int,
        vj::Int,
        Nci::Vector{Int},
        S::Int
    )

    @fastmath @inbounds begin
        r1 = 0.0   
        r2 = 0.0          
        for s = 0:2^S-1
            dig = digits(s, base = 2, pad = S)
            count = 0
            rr = 1.0
            if iseven(sum(dig))
                for nn in Nci
                    if nn != vj
                        count += 1
                        rr *= q[ci,nn,dig[count]+1]
                    end
                end
                r1 += rr
            else
                for nn in Nci
                    if nn != vj
                        count += 1
                        rr *= q[ci,nn,dig[count]+1]
                    end               
                end
                r2 += rr
            end
        end
        return r1, r2
    end
end

