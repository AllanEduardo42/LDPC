################################################################################
# Allan Eduardo Feitosa
# 26 set 2024
# Horizontal update of the LLR
# There are 4 different methods for SPA: Mckay, tanh, fast tanh and table


####################### SPA USING FAST HYPERBOLIC TANGENT ######################
function 
    calc_Lr(
        Nci::Vector{Int},
        ci::Int,
        vj::Int,    
        Lq::Matrix{Float64}
    )
    @fastmath @inbounds begin
        pLq = 1.0
        for vb in Nci
            if vb ≠ vj
                lq = Lq[ci,vb]
                if lq == 0.0
                    return 0.0
                else
                    pLq *= lq
                end
            end
        end
        if abs(pLq) < 1.0
            return 2*atanh(pLq)
        elseif signbit(pLq)
            return MINLR
        else
            return MAXLR
        end
    end
end

################ ALTERNATIVE TO HYPERBOLIC TANGENT USING TABLE ################
function 
    calc_ABCD!(
        Lq::Matrix{Float64},
        ci::Int,
        Nci::Vector{Int},
        signs::Vector{Bool},
        phi::Vector{Float64}
    )

    sLr = 0.0
    s = false
    @fastmath @inbounds for vj in Nci
        lq = Lq[ci,vj]
        sig = signbit(lq)
        s ⊻= sig
        signs[vj] = sig
        sLr += ϕ(abs(lq),phi) 
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
        lq::Float64,
        signs::Vector{Bool},
        phi::Vector{Float64}
    )

    @fastmath @inbounds begin
        x = abs(sLr - lq)
        y = signs[vj] ⊻ s
        return (1 - 2*y)*ϕ(x,phi)
    end
end

################################### MIN SUM ###################################
function 
    calc_ABCD!(
        Lq::Matrix{Float64},
        ci::Int,
        Nci::Vector{Int},
        signs::Vector{Bool},
        ::Nothing
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
        ::Float64,
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

