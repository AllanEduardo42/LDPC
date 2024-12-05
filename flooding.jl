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

        update_Lr!(Lr,Lq,m,cn2vn,Lrn,signs,phi)

    end

    # Lq update
    @inbounds for n in eachindex(vn2cn)

        _, bitvector[n] = update_Lq!(Lq,Lr,Lf[n],n,vn2cn)
        
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
        Lrn::Nothing,
        signs::Nothing,
        phi::Nothing
    )

    # horizontal update 
    @fastmath @inbounds δq = q[:,:,1]-q[:,:,2]  

    for m in eachindex(cn2vn)

        update_Lr!(r,δq,m,cn2vn)

    end

    # vertical update

    Ld = zeros(2)
    
    for n in eachindex(vn2cn)
        
        @inbounds Ld = f[n,:]
        @inbounds bitvector[n] = update_Lq!(q,r,Ld,n,vn2cn)

    end
end