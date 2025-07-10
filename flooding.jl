################################################################################
# Allan Eduardo Feitosa
# 27 set 2024
# Flooding Sum-Product Algorithm

include("update_Lr.jl")
include("update_Lq.jl")

function
    flooding!(
        bitvector::Vector{Bool},
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},
        Lf::Vector{Float64},
        Nc::Vector{Vector{Int}},
        Nv::Vector{Vector{Int}},
        signs::Union{Vector{Bool},Nothing},
        phi::Union{Vector{Float64},Nothing},
        raw::Bool
    )

    @inbounds @fastmath if raw
        # Lr update
        for ci in eachindex(Nc)
            Nci = Nc[ci]
            for vj in Nci
                Lr[ci,vj] = calc_Lr(Nci,ci,vj,Lq)
            end
        end

        # Lq update
        for vj in eachindex(Nv)
            Nvj = Nv[vj]
            Ld = calc_Ld(vj,Nvj,Lf,Lr)
            bitvector[vj] = signbit(Ld)
            for ci in Nvj
                li = LinearIndices(Lq)[ci,vj]
                Lq[li] = tanh(0.5*(Ld - Lr[li]))
            end
        end 
    else
        # Lr update
        for ci in eachindex(Nc)
            Nci = Nc[ci]
            A, B, C, D = calc_ABCD!(Lq,ci,Nci,signs,phi)
            for vj in Nci
                li = LinearIndices(Lr)[ci,vj]
                Lr[li] = calc_Lr(A,B,C,D,vj,Lq[li],signs,phi)
            end
        end

        # Lq update
        for vj in eachindex(Nv)
            Nvj = Nv[vj]
            Ld = calc_Ld(vj,Nvj,Lf,Lr)
            bitvector[vj] = signbit(Ld)
            for ci in Nvj
                li = LinearIndices(Lq)[ci,vj]
                Lq[li] = tanhLq(Ld - Lr[li],signs)
            end        
        end
    end
end

### if mode == "MKAY"

function
    flooding!(
        bitvector::Vector{Bool},
        q::Array{Float64,3},
        r::Array{Float64,3},
        f::Matrix{Float64},
        Nc::Vector{Vector{Int}},
        Nv::Vector{Vector{Int}},
        δq::Vector{Float64},
        ::Nothing,
        ::Nothing
    )

    @fastmath @inbounds begin

        # horizontal update

        for ci in eachindex(Nc)
            Nci = Nc[ci]
            # S = length(Nci)-1
            for vj in Nci
                δq[vj] = q[ci,vj,1] - q[ci,vj,2]
            end
            for vj in Nci
                δr = calc_δr(Nci,vj,δq)
                r[ci,vj,1] = 0.5*(1+δr)
                r[ci,vj,2] = 0.5*(1-δr)
                # r[ci,vj,1], r[ci,vj,2] = calc_r(q,ci,vj,Nci,S)
            end

        end
        
        # vertical update        
        for vj in eachindex(Nv)
            Nvj = Nv[vj]    
            for ci in Nvj
                Ld1, Ld2 = calc_Ld(r,f,ci,vj,Nvj)
                a = Ld1 + Ld2
                q[ci,vj,1] = Ld1/a
                q[ci,vj,2] = Ld2/a
            end
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