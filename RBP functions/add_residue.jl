################################################################################
# Allan Eduardo Feitosa
# 2 Mar 2025
# Add the residue

# List-RBP
function add_residue!(
        inlist::Union{Matrix{Bool},Nothing},
        residues::Vector{Float64},
        coords::Matrix{Int},
        residue::Float64,
        li::Int,
        ci::Int,
        vj::Int,
        listsize::Int
    )

    @fastmath @inbounds if residue > residues[listsize]

        if listsize > 1
            if residue ≥ residues[1]
                i = 1
            else
                d = listsize >> 1
                i = d
                while d > 1
                    d >>= 1
                    if residue ≥ residues[i]
                        i -= d
                    else
                        i += d
                    end
                end
                if residue < residues[i]
                    i += 1
                end
            end

            update_inlist!(inlist,coords,li)

            for j=listsize:-1:i+1
                residues[j] = residues[j-1]
                coords[1,j] = coords[1,j-1]
                coords[2,j] = coords[2,j-1]
                coords[3,j] = coords[3,j-1]
            end
        else
            i = 1
        end

        coords[1,i] = ci
        coords[2,i] = vj
        coords[3,i] = li
        residues[i] = residue        
    end
end

# RBP
function add_residue!(
    inlist::Matrix{Int},
    residues::Vector{Float64},
    ::Matrix{Int},
    residue::Float64,
    li::Int,
    ::Int,
    ::Int,
    ::Int
)

    @inbounds residues[inlist[li]] = residue

end

# auxiliary function
function update_inlist!(
    inlist::Matrix{Bool},
    coords::Matrix{Int},
    li::Int
)

    @inbounds begin
        last = coords[3,end-1]
        if last ≠ 0
            inlist[last] = false
        end
        inlist[li] = true

    end

end

function update_inlist!(
    ::Nothing,
    ::Matrix{Int},
    ::Int
)

end
