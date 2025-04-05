################################################################################
# Allan Eduardo Feitosa
# 27 set 2024
# Flooding Sum-Product Algorithm

include("update_Lr.jl")
include("update_Lq.jl")

function
    flooding!(
        bitvector::Vector{Bool},
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},
        Lf::Vector{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        Lrn::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{<:AbstractFloat},Nothing}
    )

    # Lr update
    @inbounds for m in eachindex(cn2vn)

        update_Lr!(Lr,Lq,m,cn2vn[m],Lrn,signs,phi)

    end

    # Lq update
    @inbounds for n in eachindex(vn2cn)

        bitvector[n] = update_Lq!(Lq,Lr,Lf[n],n,vn2cn[n],Lrn)
        
    end
end

### if mode == "MKAY"

function
    flooding!(
        bitvector::Vector{Bool},
        q::Array{<:AbstractFloat,3},
        r::Array{<:AbstractFloat,3},
        f::Matrix{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        δq::Vector{<:AbstractFloat},
        ::Nothing,
        ::Nothing
    )

    @inbounds begin

        # horizontal update

        for m in eachindex(cn2vn)

            vns = cn2vn[m]
            # update_Lr!(r,q,m,vns)
            for n in vns
                δq[n] = q[m,n,1] - q[m,n,2]
            end
            update_Lr!(r,δq,m,vns)

        end
        # vertical update
        
        for n in eachindex(vn2cn)    
            cns = vn2cn[n]
            update_Lq!(q,r,f[n,:],n,cns)
        end

        for n in eachindex(vn2cn) 
            d0 = f[n,1]
            d1 = f[n,2]
            for m in vn2cn[n]
                d0 *= r[m,n,1]
                d1 *= r[m,n,2]
            end
            bitvector[n] = d1 > d0
        end
    end
end