################################################################################
# Allan Eduardo Feitosa
# 12 Jan 2026
# Reliability Profile (RP) Based Low-Complexity Dynamic Schedule
# R. Yuan, T. Xie, and Z. Wang, ‘‘A reliability profile based low-complexity dynamic schedule LDPC decoding,’’ IEEE Access, vol. 10, pp. 3390–3399, 2022.

function 
    RPD!(
        bitvector::Vector{Bool},
        Lq::Matrix{Float64},
        Lr::Matrix{Float64},
        Lf::Vector{Float64},
        Nc::Vector{Vector{Int}},
        Nv::Vector{Vector{Int}},
        phi::Union{Vector{Float64},Nothing},
        msum_factor::Union{Float64,Nothing},
        oldLr::Matrix{Float64},
        syndrome::Vector{Bool},
        f::Vector{Int}
    )

    calc_syndrome!(syndrome,bitvector,Nc)

    calc_realibity!(f,Nv,syndrome)
        
    for n in eachindex(Nv)
        max_reli = 0
        vjmax = 0
        for vj in eachindex(Nv)
            reli = f[vj]
            if reli ≥ max_reli
                max_reli = reli
                vjmax = vj
            end
        end
        f[vjmax] = -1

        Nvjmax = Nv[vjmax]

        for ci in Nvjmax
            li = LinearIndices(Lq)[ci,vjmax]
            # V2C update
            Ld = calc_Ld(vjmax,Nvjmax,Lf,Lr)
            Lq[li] = Ld
            for ca in Nvjmax
                if ca < ci
                    Lq[li] -= Lr[ca,vjmax]
                elseif ca > ci
                    Lq[li] -= oldLr[ca,vjmax]
                end
            end
            # C2V update
            prod = 1.0
            Nci = Nc[ci]
            for vb in Nci
                if vb ≠ vjmax
                    prod *= tanh(Lq[ci,vb])
                end
            end
            Lr[li] = prod
            # APP
            for ca in Nvjmax
                if ca ≠ ci
                    Ld += Lr[ca,vjmax]
                end
            end
            bitvector[vjmax] = signbit(Ld)
        end
    end
end

function
    calc_realibity!(
        f::Vector{Int},
        Nv::Vector{Vector{Int}},
        syndrome::Vector{Bool}
    )

    for vj in eachindex(Nv)
        f[vj] = 0
        for ci in Nv[vj]
            f[vj] += syndrome[ci]
        end
    end
end