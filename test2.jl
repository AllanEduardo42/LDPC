using Plots
using Statistics

include("PEG.jl")
include("auxiliary_functions.jl")

function
    find_maxresidue_coords_1!(
        maxcoords::Vector{<:Integer},
        Residues::Matrix{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer}
    )

    maxresidue = 0
    for m in eachindex(cn2vn)
        for n in cn2vn[m]
            maxresidue = n
        end
    end

    return maxresidue
end

function
    find_maxresidue_coords_2!(
        maxcoords::Vector{<:Integer},
        Residues::Matrix{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer}
    )

    maxresidue = 0
    for m in eachindex(cn2vn)
        for n in @inbounds cn2vn[m]
            maxresidue = n
        end
    end

    return maxresidue
end

M,N = 256,512
D = rand([2,3,4],N)
H,_ = PEG!(D,M)

vn2cn = make_vn2cn_list(H)
cn2vn = make_cn2vn_list(H)

REALS = 1_000_000
t1 = zeros(REALS)
t2 = zeros(REALS)
Residues = zeros(M,N)
maxcoords = Vector{Int}(undef,2)

for i in 1:REALS
    Residues[H] = randn(sum(H))
    t = @timed maxresidue = find_maxresidue_coords_1!(maxcoords,Residues,cn2vn)
    t1[i] = t.time
end

for i in 1:REALS
    Residues[H] = randn(sum(H))
    t = @timed maxresidue = find_maxresidue_coords_1!(maxcoords,Residues,cn2vn)
    t2[i] = t.time
end

plotlyjs()

# histogram(
#     [t1 t3],
#     layout=grid(2,1),
#     barwidth=5e-7,
#     xlim=(0.0,0.00002),
#     ylim=(0, 10000)
# )

display((mean(t1[2:end]),mean(t2[2:end])))
display((std(t1[2:end]),std(t2[2:end])))