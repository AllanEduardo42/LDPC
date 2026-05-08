# simcore
# MKAY compatibility
# if test && bptype == "MKAY"
#     retr = Matrix{Float64}(undef,M,N)
#     retq = Matrix{Float64}(undef,M,N)
#     for ci in eachindex(Nc)
#         for vj in Nc[ci]
#             retr[ci,vj] = log.(C2V[ci,vj,1]) - log.(C2V[ci,vj,2])
#             retq[ci,vj] = log.(V2C[ci,vj,1]) - log.(V2C[ci,vj,2])
#         end
#     end
#     C2V = retr
#     V2C = retq
# end

function resetmatrix!(
    X::Array{<:Real,3},
    Nv::Vector{Vector{Int}},
    value::Real
)    
    
    @inbounds for vj in eachindex(Nv)
        for ci in Nv[vj]
            X[ci,vj,1] = value
            X[ci,vj,2] = value
        end
    end
end

# MKAY
function calc_prior_LLRs!(
    f::Matrix{Float64},
    twoZc::Int,
    signal::Vector{Float64},
    variance::Float64
)

    @fastmath @inbounds begin

        k = 1/sqrt(2*π*variance)

        for i in eachindex(signal)
            f[twoZc+i,1] = k*exp(-(signal[i]+1)^2/(2*variance))
            f[twoZc+i,2] = k*exp(-(signal[i]-1)^2/(2*variance))
        end

        normalize!(f)
    end

end

function normalize!(f::AbstractMatrix{Float64})
    
    @inbounds @fastmath for n in axes(f,1)
        α = f[n,1] + f[n,2]
        f[n,1] = f[n,1]/α
        f[n,2] = f[n,2]/α
    end

end

########################### SPA USING MKAY's METHOD ############################

function calc_post_LLR(
    r::Array{Float64,3},
    f::Matrix{Float64},
    ci::Int,
    vj::Int,
    Nvj::Vector{Int}
)

    @fastmath @inbounds begin
        Ld1 = f[vj,1]
        Ld2 = f[vj,2]
        for ca in Nvj
            if ca ≠ ci
                Ld1 *= r[ca,vj,1]
                Ld2 *= r[ca,vj,2]
            end
        end
        return Ld1, Ld2
    end
end

function init_V2C!(
    V2C::Array{Float64,3},
    prior_LLRs::Matrix{Float64},
    Nv::Vector{Vector{Int}}
)
   
    @inbounds for vj in eachindex(Nv)
        for ci in Nv[vj]
            V2C[ci,vj,1] = prior_LLRs[vj,1]
            V2C[ci,vj,2] = prior_LLRs[vj,2]
        end
    end

end

### if mode == "MKAY"

function flooding!(
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

########################### SPA USING MKAY's METHOD ############################
function calc_δr(
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
function calc_r(
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