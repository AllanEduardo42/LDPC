################################################################################
# Allan Eduardo Feitosa
# 17 set 2024
# LBP Sum-Product Algorithm

include("update_Lq.jl")
include("update_Lr.jl")
include("calc_syndrome.jl")

function
    LBP!(
        bitvector::Vector{Bool},
        Lq::Matrix{<:AbstractFloat},
        Lr::Matrix{<:AbstractFloat},        
        Lf::Vector{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        Lrn::Vector{<:AbstractFloat},
        syndrome::Vector{Bool},
        Ldn::Vector{<:AbstractFloat},
        visited_vns::Vector{Bool},
        ilbp::Bool
    )

    for m in eachindex(cn2vn)
        # Lq updates       
        @fastmath @inbounds for n in cn2vn[m] # for every n in Neighborhood(m)
            _,_ = update_Lq!(Lq,Lr,Lf[n],n,vn2cn,Lrn)
        end
        # Lr updates
        update_Lr!(Lr,Lq,m,cn2vn,Lrn,nothing,nothing)
    end

    @fastmath @inbounds for n in eachindex(vn2cn)
        m = vn2cn[n][1]
        bitvector[n] = signbit(Lq[n,m] + Lr[m,n])
    end
end