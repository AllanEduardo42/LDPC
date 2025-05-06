################################################################################
# Allan Eduardo Feitosa
# 26 set 2024
# Horizontal update of the LLR
# There are 4 different methods for SPA: Mckay, tanh, alternative and table

include("minsum.jl")

####################### SPA USING FAST HYPERBOLIC TANGENT ######################
function
    update_Lr!(
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        ci::Integer,
        Nci::Vector{<:Integer},
        Lrj::Vector{<:AbstractFloat},
        ::Nothing,
        ::Nothing     
    )

    @fastmath @inbounds begin
        pLr = 1.0
        countzeros = 0
        vj0 = 0
        for vj in Nci
            x = Lq[ci,vj]
            if x == 0.0 # Lr[ci,vj] = 0 for vj ≠ vj0
                countzeros += 1
                vj0 = vj
                Lrj[vj] = 1.0 # s.t. Lr[ci,vj0] = 2*atanh(pLr)
                if countzeros > 1 # Lr[ci,vj] = 0 ∀n
                    break
                end
            else
                Lrj[vj] = tanh(0.5*x)
            end
            pLr *= Lrj[vj]
        end
        if countzeros == 0
            for vj in Nci
                x = pLr/Lrj[vj]
                if abs(x) < 1 # controls divergent values of Lr
                    Lr[ci,vj] = 2*atanh(x)
                elseif x > 0
                    Lr[ci,vj] = INFFLOAT
                else
                    Lr[ci,vj] = NINFFLOAT
                end
            end
        else
            for vj in Nci
                Lr[ci,vj] = 0.0
            end
            if countzeros == 1 # Lr[ci,vj] = 0 for vj ≠ vj0
                Lr[ci,vj0] = 2*atanh(pLr)
            end
        end
    end
end

###################### SPA USING HYPERBOLIC TANGENT NO OPT #####################

function 
    update_Lr!(
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        ci::Integer,
        Nci::Vector{<:Integer},
        ::Nothing,
        ::Nothing,
        ::Nothing
    )
    
    @inbounds for vj in Nci
        Lr[ci,vj] = calc_Lr(Nci,ci,vj,Lq)
    end
end

function calc_Lr(
    Nci::Vector{<:Integer},
    ci::Integer,
    vj::Integer,    
    Lq::Matrix{<:AbstractFloat}
)

    @fastmath @inbounds begin
        pLr = 1.0
        for vb in Nci
            if vb ≠ vj
                pLr *= tanh(0.5*Lq[ci,vb])
            end
        end
        return 2*atanh(pLr)
    end
end


################# ALTERNATIVE TO HYPERBOLIC TANGENT AND TABLE ##################

function ϕ(Lq::AbstractFloat, ::Nothing)    
    @fastmath log(1 + 2/(exp(Lq)-1))
end

function ϕ(Lq::AbstractFloat, phi::Vector{<:AbstractFloat})    
    @inbounds phi[get_index(Lq)]
end

function
    phi_sign!(
        Lq::AbstractFloat,
        s::Bool,
        phi::Union{Vector{<:AbstractFloat},Nothing}
    )

    sig = signbit(Lq)
    @fastmath ab = abs(Lq)
    return ϕ(ab,phi), sig, s ⊻ sig

end

function 
    update_Lr!(
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        ci::Integer,
        Nci::Vector{<:Integer},
        Lrj::Vector{<:AbstractFloat},
        signs::Vector{Bool},
        phi::Union{Vector{<:AbstractFloat},Nothing}
    )

    @fastmath @inbounds begin
        sLr = 0.0
        s = false
        for vj in Nci
            Lrj[vj], signs[vj], s = phi_sign!(Lq[ci,vj],s,phi)
            sLr += Lrj[vj] 
        end
        for vj in Nci
            x = abs(sLr - Lrj[vj])
            if x > 0 # (Inf restriction)
                Lr[ci,vj] = (1 - 2*(signs[vj] ⊻ s))*ϕ(x,phi)
            end
        end
    end   
end

################################### MIN SUM ###################################


function
    update_Lr!(
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        ci::Integer,
        Nci::Vector{<:Integer},
        ::Nothing,
        signs::Vector{Bool},
        ::Nothing       
    )

    minsum!(Lq,Lr,signs,ci,Nci)

end


########################### SPA USING MKAY's METHOD ############################
function
    update_Lr!(
        r::Array{<:AbstractFloat,3},
        δq::Vector{<:AbstractFloat},
        ci::Integer,
        Nci::Vector{<:Integer}    
    )
    
    @inbounds for vj in Nci
        δr = 1.0
        for vb in Nci
            if vb ≠ vj
                δr *= δq[vb]
            end
        end
        r[ci,vj,1] = 0.5*(1+δr)
        r[ci,vj,2] = 0.5*(1-δr)
    end 

end

# pre-historic method
function
    update_Lr!(
        r::Array{<:AbstractFloat,3},
        q::Array{<:AbstractFloat,3},
        ci::Integer,
        Nci::Vector{<:Integer}    
    )

    S = length(Nci)-1
    @inbounds for vj in Nci
        r[ci,vj,1] = 0.0   
        r[ci,vj,2] = 0.0          
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
                r[ci,vj,1] += rr
            else
                for nn in Nci
                    if nn != vj
                        count += 1
                        rr *= q[ci,nn,dig[count]+1]
                    end               
                end
                r[ci,vj,2] += rr
            end
        end
    end
end

