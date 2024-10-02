function
    f1(
        maxcoords::Vector{<:Integer},
        Residues::Matrix{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        reals::Integer
    )
    # residue = 0.0
    # comp = false
    for r in 1:reals
        maxresidue = 0
        for m in eachindex(cn2vn)
            for n in cn2vn[m]
                @inbounds residue = Residues[m,n]
                @fastmath comp = residue > maxresidue
                if comp
                    maxresidue = residue
                    @inbounds maxcoords[1] = m
                    @inbounds maxcoords[2] = n
                end
            end
        end
    end
end

REALS=100000
vn2cn = make_vn2cn_list(H)
cn2vn = make_cn2vn_list(H)

Residues = zeros(M,N)
maxcoords = Vector{Int}(undef,2)


Residues[H] = randn(sum(H))
