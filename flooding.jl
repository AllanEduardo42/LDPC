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
        Nc::Vector{Vector{T}} where {T<:Integer},
        Nv::Vector{Vector{T}} where {T<:Integer},
        Lrj::Union{Vector{<:AbstractFloat},Nothing},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{<:AbstractFloat},Nothing}
    )

    # Lr update
    @inbounds for ci in eachindex(Nc)

        update_Lr!(Lr,Lq,ci,Nc[ci],Lrj,signs,phi)

    end

    # Lq update
    @inbounds for vj in eachindex(Nv)

        bitvector[vj] = update_Lq!(Lq,Lr,Lf,vj,Nv[vj],Lrj)
        
    end
end

### if mode == "MKAY"

function
    flooding!(
        bitvector::Vector{Bool},
        q::Array{<:AbstractFloat,3},
        r::Array{<:AbstractFloat,3},
        f::Matrix{<:AbstractFloat},
        Nc::Vector{Vector{T}} where {T<:Integer},
        Nv::Vector{Vector{T}} where {T<:Integer},
        δq::Vector{<:AbstractFloat},
        ::Nothing,
        ::Nothing
    )

    @inbounds begin

        # horizontal update

        for ci in eachindex(Nc)

            vns = Nc[ci]
            # update_Lr!(r,q,ci,vns)
            for vj in vns
                δq[vj] = q[ci,vj,1] - q[ci,vj,2]
            end
            update_Lr!(r,δq,ci,vns)

        end
        
        # vertical update
        
        for vj in eachindex(Nv)    
            update_Lq!(q,r,f[vj,:],vj,Nv[vj])
        end

        for vj in eachindex(Nv) 
            d0 = f[vj,1]
            d1 = f[vj,2]
            for ci in Nv[vj]
                d0 *= r[ci,vj,1]
                d1 *= r[ci,vj,2]
            end
            bitvector[vj] = d1 > d0
        end
    end
end