using Plots
# using FastRunningMedian

include("auxiliary_functions.jl")
include("PEG.jl")

println(Threads.nthreads())

function
    find_maxresidue_coords_multi!(
        Residues::Matrix{<:AbstractFloat},
        maxcoordsth::Matrix{<:Integer},
        chunks,
        offsets
    )
    tasks = map(chunks,offsets) do chunk,offset
        Threads.@spawn find_maxresidue_coords!(Residues,view(maxcoordsth,:,Threads.threadid()),chunk,offset)
    end
    maxresidues = fetch.(tasks)
    a,b = findmax(maxresidues)
    return a, maxcoordsth[:,b]
end

function
    find_maxresidue_coords!(
        Residues::Matrix{<:AbstractFloat},
        maxcoords::AbstractVector{<:Integer},
        cn2vn::AbstractVector{Vector{T}} where {T<:Integer},
        offset::Integer,
    )

    maxresidue = 0.0
    for m in eachindex(cn2vn)
        for n in cn2vn[m]
            @inbounds residue = Residues[m+offset,n]
            if @fastmath residue > maxresidue
                maxresidue = residue
                maxcoords[1] = m+offset
                maxcoords[2] = n
            end
        end
    end

    return maxresidue
end

M = 512
N = 1024
SEED_GRAPH::Int64 = 5714
rng_graph = Xoshiro(SEED_GRAPH)
D = rand(rng_graph,[2,3,4],N)
if !@isdefined(H) 
    H,_ = PEG(D,M)
end
cn2vn  = make_cn2vn_list(H)
Residues = 0.0*H
Residues[H] = randn(sum(H))
maxcoords = [0,0]

num = length(cn2vn) ÷ Threads.nthreads()
chunks = Iterators.partition(cn2vn,num)
offsets = num*(0:Threads.nthreads()-1)
maxcoordsth = zeros(Int,2,Threads.nthreads())

@time a = find_maxresidue_coords!(Residues,maxcoords,cn2vn,0)
@time b = find_maxresidue_coords_multi!(Residues,maxcoordsth,chunks,offsets)

display((a,maxcoords))
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