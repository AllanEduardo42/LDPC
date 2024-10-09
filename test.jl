Residues = zeros(M,N)

num_edges = sum(H)

Residues[H] = randn(num_edges)

Ii = sortperm(Residues[H],rev=true)

Ii_inv = sortperm(Ii)

r = Residues[H][Ii]

Addrs_inv = zeros(Int,M,N)
Addrs = zeros(Int,2,num_edges)

vn2cn = make_vn2cn_list(H)

e = 0
for n in eachindex(vn2cn)
    for m in vn2cn[n]
        global e +=1
        Addrs_inv[m,n] = e
        Addrs[1,e] = m
        Addrs[2,e] = n
    end
end

REALS = 10000

MIN = true

for real in 1:REALS
    local position = collect(1:num_edges)
    n = rand(1:N)
    m = rand(vn2cn[n])
    Residues[m,n] = randn()
    local idx = Ii_inv[Addrs_inv[m,n]]
    r[idx] = Residues[m,n]
    for e in eachindex(r)
        if r[idx] > r[e]
            if idx > e            
                position[e:idx-1] .+= 1
                position[idx] = e
            else
                position[idx] = e-1
                position[idx+1:e-1] .-= 1              
            end
            local v = sortperm(position)
            global Ii = Ii[v]
            global Ii_inv = sortperm(Ii)
            global r = r[v]
            global MIN = false
            if r != sort(Residues[H],rev=true)
                println("break 1 at real = ", real)
                display((r .!= sort(Residues[H],rev=true))')
                break
            end
            break
        end
    end
    if MIN
        position[idx] = num_edges
        position[idx+1:end] .-= 1
        local v = sortperm(position)
        global Ii = Ii[v]
        global Ii_inv = sortperm(Ii)
        global r = r[v]
        if r != sort(Residues[H],rev=true)
            println("break 2 at real = ", real)
            display((r .!= sort(Residues[H],rev=true))')
            break
        end
    end 
    # if Ii != sortperm(Residues[H],rev=true)
    #     println("break 1 at real = ", real)
    #     break
    # end
end
