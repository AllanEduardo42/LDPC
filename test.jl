using Plots
using FastRunningMedian

include("auxiliary_functions.jl")
include("PEG.jl")

Threads.nthreads()

function sum_single(a)
    s = 0
    for r in a
        s += r
    end
    s
end

function sum_multi_good(a)
    chunks = Iterators.partition(a, length(a) ÷ Threads.nthreads())
    tasks = map(chunks) do chunk
        Threads.@spawn sum_single(chunk)
    end
    chunk_sums = fetch.(tasks)
    return sum_single(chunk_sums)
end

function
    find_maxresidue_coords_multi!(
        Residues::Matrix{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer}
    )
    num = length(cn2vn) ÷ Threads.nthreads()
    chunks = Iterators.partition(cn2vn,num)
    offsets = num*(0:Threads.nthreads()-1)
    tasks = map(chunks,offsets) do chunk,offset
        Threads.@spawn find_maxresidue_coords!(Residues,chunk,offset)
    end
    chunk_ = fetch.(tasks)
    return maximum(chunk_)
end

function
    find_maxresidue_coords!(
        Residues::Matrix{<:AbstractFloat},
        cn2vn::AbstractVector{Vector{T}} where {T<:Integer},
        offset::Integer,
    )

    maxresidue = 0.0
    M = 0
    N = 0
    for m in eachindex(cn2vn)
        for n in cn2vn[m]
            @inbounds residue = Residues[m+offset,n]
            if @fastmath residue > maxresidue
                maxresidue = residue
                M = m
                N = n
            end
        end
    end

    return maxresidue, M+offset, N
end

M = 512
N = 1024
SEED_GRAPH::Int64 = 5714
rng_graph = Xoshiro(SEED_GRAPH)
D = rand(rng_graph,[2,3,4],N)
H,_ = PEG(D,M)
cn2vn  = make_cn2vn_list(H)
Residues = 0.0*H
Residues[H] = randn(sum(H))

@time a = find_maxresidue_coords!(Residues,cn2vn,0)
@time b = find_maxresidue_coords_multi!(Residues,cn2vn)

display(a)
display(b)

# R = 10000
# t1 = zeros(length(10:10:R))
# t2 = zeros(length(10:10:R))

# j = 0
# for r = 10:10:R
#     global j += 1
#     x = randn(r)
#     t = @timed sum_single(x)
#     t1[j] = t.time
# end

# j = 0
# for r = 10:10:R
#     global j += 1
#     x = randn(r)
#     t = @timed sum_multi_good(x)
#     t2[j] = t.time

# end

# plotlyjs()

# plot(20:10:R,running_median(t1[2:end],3),label="no threading")
# plot!(20:10:R,running_median(t2[2:end],3),label="with threading")