################################################################################
# Allan Eduardo Feitosa
# 26 set 2024
# Horizontal update of the LLR
# There are 4 different methods for SPA: Mckay, tanh, alternative and table

include("minsum.jl")

########################### SPA USING MKAY's METHOD ############################
function
    update_Lr!(
        r::Array{<:AbstractFloat,3},
        δq::Matrix{<:AbstractFloat},
        m::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer},     
    )
    
    for n in cn2vn[m]
        δr = 1
        for n2 in cn2vn[m]
            if n2 ≠ n
                @fastmath @inbounds δr *= δq[n2,m]
            end
        end
        @fastmath @inbounds r[m,n,1] = 0.5*(1+δr)
        @fastmath @inbounds r[m,n,2] = 0.5*(1-δr)
    end
   
        

end

######################### SPA USING HYPERBOLIC TANGENT #########################
function
    update_Lr!(
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        m::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        Lrn::Vector{<:AbstractFloat},
        ::Nothing,
        ::Nothing        
    )
    pLr = 1.0
    for n in cn2vn[m]
        @fastmath @inbounds Lrn[n] = tanh(0.5*Lq[n,m])
        @fastmath @inbounds pLr *= Lrn[n]
    end
    for n in cn2vn[m]
        @fastmath @inbounds x = pLr/Lrn[n]
        if @fastmath abs(x) < 1 # controls divergent values of Lr
            @fastmath @inbounds Lr[m,n] = 2*atanh(x)
        end
    end
end

###################### SPA USING HYPERBOLIC TANGENT VER 2 ######################

function 
    update_Lr!(
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        m::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        ::Nothing,
        ::Nothing,
        ::Nothing   
    )

    for n in cn2vn[m]
        @inbounds Lr[m,n] = 1.0
        for n2 in cn2vn[m]
            if n2 ≠ n
                @fastmath @inbounds Lr[m,n] *= tanh(0.5*Lq[n2,m])
            end
        end
        if @fastmath @inbounds abs(Lr[m,n]) < 1
            @fastmath @inbounds Lr[m,n] = 2*atanh(Lr[m,n])
        else
            @fastmath @inbounds Lr[m,n] = Lr[m,n]*INF
        end
    end
end


################# ALTERNATIVE TO HYPERBOLIC TANGENT AND TABLE ##################

function ϕ(Lq::AbstractFloat, ::Nothing)    
    @fastmath log(1 + 2/(exp(Lq)-1))
end

function ϕ(Lq::AbstractFloat, phi::Vector{<:AbstractFloat})    
    phi[get_index(Lq)]
end

function
    phi_sign!(
        Lq::AbstractFloat,
        s::Bool,
        phi::Union{Vector{<:AbstractFloat},Nothing}
    )

    _sign = signbit(Lq)
    return ϕ(abs(Lq),phi), _sign, s ⊻ _sign

end

function 
    update_Lr!(
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        m::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        Lrn::Vector{<:AbstractFloat},
        signs::Vector{Bool},
        phi::Union{Vector{<:AbstractFloat},Nothing}
    )

    sLr = 0.0
    s = false
    for n in cn2vn[m]
        @inbounds Lrn[n], signs[n], s = phi_sign!(Lq[n,m],s,phi)
        @fastmath @inbounds sLr += Lrn[n] 
    end
    for n in cn2vn[m]
        @fastmath @inbounds x = abs(sLr - Lrn[n])
        if @fastmath x > 0 # (Inf restriction)
            @inbounds Lr[m,n] = (1 - 2*(signs[n] ⊻ s))*ϕ(x,phi)
        end
    end    
end

################################### MIN SUM ###################################


function
    update_Lr!(
        Lr::Matrix{<:AbstractFloat},
        Lq::Matrix{<:AbstractFloat},
        m::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        ::Nothing,
        signs::Vector{Bool},
        ::Nothing        
    )

    minsum!(Lq,Lr,signs,m,cn2vn)

end
