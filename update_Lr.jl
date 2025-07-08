################################################################################
# Allan Eduardo Feitosa
# 26 set 2024
# Horizontal update of the LLR
# There are 4 different methods for SPA: Mckay, tanh, fast tanh and table


####################### SPA USING FAST HYPERBOLIC TANGENT ######################
function 
    calc_ABCD!(
        aux::Vector{Float64},
        ::Nothing,
        ::Nothing,
        Lq::Matrix{Float64},
        ci::Int,
        Nci::Vector{Int},
    )

    pLr = 1.0
    count_zeros = 0
    vj_notzero = 0
    @fastmath @inbounds for vj in Nci
        lq = Lq[ci,vj]
        if lq == 0.0 # Lr[ci,vj] = 0 for vj ≠ vj0
            if count_zeros == 1 # Lr[ci,vj0] = 0
                count_zeros = 2
                break
            end
            count_zeros = 1
            vj_notzero = vj
        else
            aux[vj] = lq
            pLr *= lq
        end
    end
    
    return pLr, count_zeros, vj_notzero, nothing
end

function 
    calc_Lr(
        pLr::Float64,           #A
        count_zeros::Int,       #B
        vj_notzero::Int,        #C
        ::Nothing,              #D
        vj::Int,
        aux::Vector{Float64},
        ::Nothing,
        ::Nothing
    )

    @fastmath @inbounds begin
        if count_zeros == 0
            x = pLr/aux[vj]
            if abs(x) < 1 # controls divergent values of Lr
                return 2*atanh(x)
            elseif x > 0
                return INFFLOAT
            else
                return NINFFLOAT
            end
        elseif count_zeros == 1 && vj_notzero == vj
            return 2*atanh(pLr)
        else
            return 0.0
        end
    end
end

###################### SPA USING HYPERBOLIC TANGENT NO OPT #####################
function 
    calc_Lr(
        Nci::Vector{Int},
        ci::Int,
        vj::Int,    
        Lq::Matrix{Float64}
    )

    @fastmath @inbounds begin
        pLr = 1.0
        for vb in Nci
            if vb ≠ vj
                lq = Lq[ci,vb]
                if lq == 0.0
                    return 0.0
                else
                    pLr *= tanh(0.5*lq)
                end
            end
        end
        return 2*atanh(pLr)
    end
end


################ ALTERNATIVE TO HYPERBOLIC TANGENT USING TABLE ################
function 
    calc_ABCD!(
        aux::Vector{Float64},
        signs::Vector{Bool},
        phi::Vector{Float64},
        Lq::Matrix{Float64},
        ci::Int,
        Nci::Vector{Int},
    )

    sLr = 0.0
    s = false
    @fastmath @inbounds for vj in Nci
        lq = Lq[ci,vj]
        sig = signbit(lq)
        s ⊻= sig
        aux[vj] = ϕ(abs(lq),phi)
        signs[vj] = sig
        sLr += aux[vj] 
    end

    return sLr, s, nothing, nothing
end

function 
    calc_Lr(
        sLr::Float64,       #A
        s::Bool,            #B          
        ::Nothing,          #C
        ::Nothing,          #D
        vj::Int,
        aux::Vector{Float64},
        signs::Vector{Bool},
        phi::Vector{Float64}
    )

    @fastmath @inbounds begin
        x = abs(sLr - aux[vj])
        y = signs[vj] ⊻ s
        return (1 - 2*y)*ϕ(x,phi)
    end
end

################################### MIN SUM ###################################
function 
    calc_ABCD!(
        ::Vector{Float64},
        signs::Vector{Bool},
        ::Nothing,
        Lq::Matrix{Float64},
        ci::Int,
        Nci::Vector{Int}
    )

    @fastmath @inbounds begin
        s = false
        minL = INFFLOAT
        minL2 = INFFLOAT
        vjmin = Nci[1]
        for vj in Nci
            lq = Lq[ci,vj]
            sig = signbit(lq)
            s ⊻= sig
            β = abs(lq)
            signs[vj] = sig
            if β < minL
                vjmin = vj
                minL, minL2 = β, minL
            elseif β < minL2
                minL2 = β
            end
        end
    end

    return minL*ALPHA, s, vjmin, minL2*ALPHA

end

function 
    calc_Lr(
        minL::Float64,
        s::Bool,
        vjmin::Int,
        minL2::Float64,
        vj::Int,
        ::Vector{Float64},
        signs::Vector{Bool},
        ::Nothing
    )

    @fastmath @inbounds begin
        if signs[vj] ⊻ s
            if vj ≠ vjmin
                return -minL
            else
                return -minL2
            end
        else
            if vj ≠ vjmin
                return minL
            else
                return minL2
            end
        end
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

