################################################################################
# Allan Eduardo Feitosa
# 26 set 2024
# Horizontal update of the LLR
# There are 4 different methods for SPA: Mckay, tanh, alternative and table

include("min_sum.jl")

########################### SPA USING MKAY's METHOD ############################
function
    update_Lr!(
        r::Array{<:AbstractFloat,3},
        δq::Matrix{<:AbstractFloat},
        m::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer},     
    )
    δr = 1
    for n in cn2vn[m]
        @inbounds δr *= δq[n,m]
    end
    for n in cn2vn[m]
        x = δr/(δq[n,m] + eps())        
        @inbounds r[m,n,1] = 0.5*(1+x)
        @inbounds r[m,n,2] = 0.5*(1-x)
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
        signs::Nothing,
        phi::Nothing        
    )
    pLr = 1.0
    for n in cn2vn[m]
        @inbounds @fastmath Lrn[n] = tanh(0.5*Lq[n,m])
        @inbounds @fastmath pLr *= Lrn[n]
    end
    for n in cn2vn[m]
        @inbounds @fastmath x = pLr/Lrn[n]
        if abs(x) < 1 # controls divergent values of Lr
            @inbounds @fastmath Lr[m,n] = 2*atanh(x)
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
        Lrn::Nothing,
        signs::Nothing,
        phi::Nothing   
    )

    for n in cn2vn[m]
        @inbounds Lr[m,n] = 1.0
        for n2 in cn2vn[m]
            if n2 ≠ n
                @inbounds @fastmath Lr[m,n] *= tanh(0.5*Lq[n2,m])
            end
        end
        @inbounds @fastmath Lr[m,n] = 2*atanh(Lr[m,n])
    end
end


################# ALTERNATIVE TO HYPERBOLIC TANGENT AND TABLE ##################

function ϕ(Lq::AbstractFloat, phi::Nothing)    
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
        @inbounds @fastmath sLr += Lrn[n] 
    end
    for n in cn2vn[m]
        x = abs(sLr - Lrn[n])
        if x > 0 # (Inf restriction)
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
        Lrn::Nothing,
        signs::Vector{Bool},
        phi::Nothing        
    )
    args = _min_sum!(Lq,signs,m,cn2vn)
    for n in cn2vn[m]
        Lr[m,n] = __min_sum!(n,signs[n],args...)
    end 

end
