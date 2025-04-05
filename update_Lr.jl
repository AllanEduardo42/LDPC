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
        m::Integer,
        vns::Vector{<:Integer},
        Lrn::Vector{<:AbstractFloat},
        ::Nothing,
        ::Nothing     
    )

    @fastmath @inbounds begin
        pLr = 1.0
        countzeros = 0
        n0 = 0
        for n in vns
            x = Lq[m,n]
            if x == 0.0 # Lr[m,n] = 0 for n ≠ n0
                countzeros += 1
                n0 = n
                Lrn[n] = 1.0 # s.t. Lr[m,n0] = 2*atanh(pLr)
                if countzeros > 1 # Lr[m,n] = 0 ∀n
                    break
                end
            else
                Lrn[n] = tanh(0.5*x)
            end
            pLr *= Lrn[n]
        end
        if countzeros == 0
            for n in vns
                x = pLr/Lrn[n]
                if abs(x) < 1 # controls divergent values of Lr
                    Lr[m,n] = 2*atanh(x)
                elseif x > 0
                    Lr[m,n] = INFFLOAT
                else
                    Lr[m,n] = NINFFLOAT
                end
            end
        else
            for n in vns
                Lr[m,n] = 0.0
            end
            if countzeros == 1 # Lr[m,n] = 0 for n ≠ n0
                Lr[m,n0] = 2*atanh(pLr)
            end
        end
    end
end

###################### SPA USING HYPERBOLIC TANGENT NO OPT #####################

function 
    update_Lr!(
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        m::Integer,
        vns::Vector{<:Integer},
        ::Nothing,
        ::Nothing,
        ::Nothing
    )
    
    @inbounds for n in vns
        pLr = 1.0
        for n2 in vns
            if n2 ≠ n
                pLr *= tanh(0.5*Lq[m,n2])
            end
        end
        Lr[m,n] = 2*atanh(pLr)
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
        m::Integer,
        vns::Vector{<:Integer},
        Lrn::Vector{<:AbstractFloat},
        signs::Vector{Bool},
        phi::Union{Vector{<:AbstractFloat},Nothing}
    )

    @fastmath @inbounds begin
        sLr = 0.0
        s = false
        for n in vns
            Lrn[n], signs[n], s = phi_sign!(Lq[m,n],s,phi)
            sLr += Lrn[n] 
        end
        for n in vns
            x = abs(sLr - Lrn[n])
            if x > 0 # (Inf restriction)
                Lr[m,n] = (1 - 2*(signs[n] ⊻ s))*ϕ(x,phi)
            end
        end
    end   
end

################################### MIN SUM ###################################


function
    update_Lr!(
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        m::Integer,
        vns::Vector{<:Integer},
        ::Nothing,
        signs::Vector{Bool},
        ::Nothing       
    )

    minsum!(Lq,Lr,signs,m,vns)

end


########################### SPA USING MKAY's METHOD ############################
function
    update_Lr!(
        r::Array{<:AbstractFloat,3},
        δq::Vector{<:AbstractFloat},
        m::Integer,
        vns::Vector{<:Integer}    
    )
    
    @inbounds for n in vns
        δr = 1.0
        for n2 in vns
            if n2 ≠ n
                δr *= δq[n2]
            end
        end
        r[m,n,1] = 0.5*(1+δr)
        r[m,n,2] = 0.5*(1-δr)
    end 

end

# pre-historic method
function
    update_Lr!(
        r::Array{<:AbstractFloat,3},
        q::Array{<:AbstractFloat,3},
        m::Integer,
        vns::Vector{<:Integer}    
    )

    S = length(vns)-1
    @inbounds for n in vns
        r[m,n,1] = 0.0   
        r[m,n,2] = 0.0          
        for s = 0:2^S-1
            dig = digits(s, base = 2, pad = S)
            count = 0
            rr = 1.0
            if iseven(sum(dig))
                for nn in vns
                    if nn != n
                        count += 1
                        rr *= q[m,nn,dig[count]+1]
                    end
                end
                r[m,n,1] += rr
            else
                for nn in vns
                    if nn != n
                        count += 1
                        rr *= q[m,nn,dig[count]+1]
                    end               
                end
                r[m,n,2] += rr
            end
        end
    end
end

