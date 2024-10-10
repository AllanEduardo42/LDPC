function
    find_maxresidue_coords_test!(
        Residues::Matrix{<:AbstractFloat},
        list::Vector{Vector{T}} where {T<:Integer},
        offset::Integer
    )

    maxresidue = 0.0
    cnmax = 0
    vnmax = 0
    for m in eachindex(list)
        for n in list[m]
            @inbounds residue = Residues[m+offset,n]
            if @fastmath residue > maxresidue
                maxresidue = residue
                cnmax = m+offset
                vnmax = n
            end
        end
    end

    return maxresidue, cnmax, vnmax
end

function 
    multi(
        Residues,
        maxresidues,
        cnmax,
        vnmax,
        lists,
        offsets,
        chunk_length
    )

    Threads.@threads for i in 1:nthreads
        maxresidues[i], cnmax[i], vnmax[i] = find_maxresidue_coords_test!(Residues,lists[i],offsets[i])
    end

    a,b = findmax(maxresidues)

    return a, cnmax[b], vnmax[b]

end

nthreads = 32
chunk_length = M÷nthreads
offsets = 0:chunk_length:M-chunk_length

Residues = 0.0*H
Residues[H] = abs.(randn(sum(H)))
cn2vn = make_cn2vn_list(H)

maxresidues = zeros(nthreads)
cnmaxs = zeros(Int,nthreads)
vnmaxs = zeros(Int,nthreads)

lists = Vector{Vector{Vector{Int}}}()
for i=1:nthreads
    push!(lists,cn2vn[offsets[i]+1:offsets[i]+chunk_length])
end

@time maxresidue, cnmax, vnmax = find_maxresidue_coords_test!(Residues,cn2vn,0)
display((maxresidue, cnmax, vnmax))
@time maxresidue, cnmax, vnmax = multi(Residues,maxresidues,cnmaxs,vnmaxs,lists,offsets,chunk_length)
display((maxresidue, cnmax, vnmax))