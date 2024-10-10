function
    RBP_test!(
        Residues::Matrix{<:AbstractFloat},
        cnmax::Integer,
        vnmax::Integer,
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        num_edges::Union{Integer,Nothing},
        Ii::Union{Vector{Int},Nothing},
        Ii_inv::Union{Vector{Int},Nothing},
        r::Union{Vector{Float64},Nothing},
        Addrs_inv::Union{Matrix{Int},Nothing},
        rng_noise
    )

    for m in vn2cn[vnmax]
        if m ≠ cnmax
            for n in cn2vn[m]
                if n ≠ vnmax
                    residue = abs(randn(rng_noise))
                    core(Residues,residue,m,n,num_edges,Ii,Ii_inv,r,Addrs_inv)
                end
            end
        end
    end

end

function 
    core(
        Residues::Matrix{<:AbstractFloat},
        residue::Float64,
        m::Integer,
        n::Integer,
        num_edges::Integer,
        Ii::Nothing,
        Ii_inv::Nothing,
        r::Nothing,
        Addrs_inv::Nothing        
    )
    return @fastmath @inbounds Residues[m,n] = residue
end

function 
    core(
        Residues::Matrix{<:AbstractFloat},
        residue::Float64,
        m::Integer,
        n::Integer,
        num_edges::Integer,
        Ii::Vector{Int},
        Ii_inv::Vector{Int},
        r::Vector{Float64},
        Addrs_inv::Matrix{Int}
    )

    MIN = true
    position = collect(1:num_edges)
    idx = Ii_inv[Addrs_inv[m,n]]   
    r[idx] = residue
    Residues[m,n] = r[idx]
    for e in 1:num_edges
        if r[idx] > r[e]
            if idx > e 
                println("UM")      
                position[e:idx-1] .+= 1
                position[idx] = e
            else
                println("DOIS")  
                position[idx] = e-1
                position[idx+1:e-1] .-= 1              
            end
            v = sortperm(position)
            Ii = Ii[v]
            print("1: ")
            diff = Ii .!= sortperm(Residues[H],rev=true)
            display(sum(diff))
            p = findfirst(x -> x == 1, diff)
            if p !== nothing 
                display(Ii[p-1:p+1])
                display(sortperm(Residues[H],rev=true)[p-1:p+1])
            end
            Ii_inv = sortperm(Ii)
            r .= r[v]
            MIN = false
            break
        end
    end
    if MIN
        position[idx] = num_edges
        position[idx+1:end] .-= 1
        v = sortperm(position)
        Ii = Ii[v]
        Ii_inv = sortperm(Ii)
        r .= r[v]
    end
end

function
    find_maxresidue_coords_test!(
        Residues::Matrix{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer}
    )

    maxresidue = 0.0
    cnmax = 0
    vnmax = 0
    for m in eachindex(cn2vn)
        for n in cn2vn[m]
            @inbounds residue = Residues[m,n]
            if @fastmath residue > maxresidue
                maxresidue = residue
                cnmax = m
                vnmax = n
            end
        end
    end

    return maxresidue, cnmax, vnmax
end

function 
    f(
        Residues::Matrix{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        cnmax::Integer,
        vnmax::Integer,
        num_edges::Integer,
        Ii::Nothing,
        Ii_inv::Nothing,
        r::Nothing,
        Addrs_inv::Nothing,
        Addrs::Nothing,
        rng_noise
    )

    @inbounds Residues[cnmax,vnmax] = 0.0
    RBP_test!(
        Residues,
        cnmax,
        vnmax,
        cn2vn,
        vn2cn,
        num_edges,
        Ii,
        Ii_inv,
        r,
        Addrs_inv,
        rng_noise
    )
    maxresidue, cnmax, vnmax = find_maxresidue_coords_test!(Residues,cn2vn) 

    return maxresidue, cnmax, vnmax

end

function 
    f(
        Residues::Matrix{<:AbstractFloat},
        cn2vn::Vector{Vector{T}} where {T<:Integer},
        vn2cn::Vector{Vector{T}} where {T<:Integer},
        cnmax::Integer,
        vnmax::Integer,
        num_edges::Integer,
        Ii::Vector{Int},
        Ii_inv::Vector{Int},
        r::Vector{Float64},
        Addrs_inv::Matrix{Int},
        Addrs::Matrix{Int},
        rng_noise
    )

    @inbounds Residues[cnmax,vnmax] = 0.0    
    idx = Ii_inv[Addrs_inv[cnmax,vnmax]]
    r[idx] = 0.0
    position = collect(1:num_edges)
    position[idx] = num_edges
    position[idx+1:end] .-= 1
    v = sortperm(position)
    Ii = Ii[v]
    print("0: ")
    display(sum(Ii .!= sortperm(Residues[H],rev=true)))
    Ii_inv = sortperm(Ii)
    r .= r[v]

    RBP_test!(
        Residues,
        cnmax,
        vnmax,
        cn2vn,
        vn2cn,
        num_edges,
        Ii,
        Ii_inv,
        r,
        Addrs_inv,
        rng_noise
    )

    return r[1], Addrs[1,Ii[1]], Addrs[2,Ii[1]]

end

function 
    loop(
        reals,
        Residues,
        cn2vn,
        vn2cn,
        cnmax,
        vnmax,
        num_edges,
        Ii,
        Ii_inv,
        r,
        Addrs_inv,
        Addrs,
        rng_noise
    )

    maxresidue = 0
    for real in 1:reals
        maxresidue, cnmax, vnmax = f(
            Residues,
            cn2vn,
            vn2cn,
            cnmax,
            vnmax,
            num_edges,
            Ii,
            Ii_inv,
            r,
            Addrs_inv,
            Addrs,
            rng_noise
        )
    end

    return maxresidue, cnmax, vnmax

end


cn2vn = make_cn2vn_list(H)
vn2cn = make_vn2cn_list(H)
Residues = zeros(M,N)
num_edges = sum(H)
Addrs_inv = zeros(Int,M,N)
Addrs = zeros(Int,2,num_edges)

e = 0
for n in eachindex(vn2cn)
    for m in vn2cn[n]
        global e +=1
        Addrs_inv[m,n] = e
        Addrs[1,e] = m
        Addrs[2,e] = n
    end
end

REALS = 1

rng_noise = Xoshiro(1428)
Residues[H] = abs.(randn(rng_noise, num_edges)) 
_, cnmax, vnmax = find_maxresidue_coords_test!(Residues,cn2vn)

@time Maxresidue, Cnmax, Vnmax = loop(
    REALS,
    Residues,
    cn2vn,
    vn2cn,
    cnmax,
    vnmax,
    num_edges,
    nothing,
    nothing,
    nothing,
    nothing,
    nothing,
    rng_noise
)
display((Maxresidue, Cnmax, Vnmax, findmax(Residues)))

rng_noise = Xoshiro(1428)
Residues[H] = abs.(randn(rng_noise, num_edges))
_, Cnmax, Vnmax = find_maxresidue_coords_test!(Residues,cn2vn) 

Ii = sortperm(Residues[H],rev=true)
Ii_inv = sortperm(Ii)    
R = Residues[H][Ii]

@time Maxresidue, Cnmax, Vnmax = loop(
    REALS,
    Residues,
    cn2vn,
    vn2cn,
    cnmax,
    vnmax,
    num_edges,
    Ii,
    Ii_inv,
    R,
    Addrs_inv,
    Addrs,
    rng_noise
)
display((Maxresidue, Cnmax, Vnmax, findmax(Residues)))
