################################################################################
# Allan Eduardo Feitosa
# 2 Mar 2025
# Add the residue

# List-RBP
function add_residue!(
        inlist::Union{Matrix{Bool},Nothing},
        residues::Vector{<:AbstractFloat},
        coords::Matrix{<:Integer},
        residue::AbstractFloat,
        li::Integer,
        m::Integer,
        n::Integer,
        listsize::Integer
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

        coords[1,i] = m
        coords[2,i] = n
        coords[3,i] = li
        residues[i] = residue        
    end
end

# RBP
function add_residue!(
    inlist::Matrix{<:Integer},
    residues::Vector{<:AbstractFloat},
    ::Matrix{<:Integer},
    residue::AbstractFloat,
    li::Integer,
    ::Integer,
    ::Integer,
    ::Integer
)

    @inbounds residues[inlist[li]] = residue

end

# auxiliary function
function update_inlist!(
    inlist::Matrix{Bool},
    coords::Matrix{<:Integer},
    li::Integer
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
    ::Matrix{<:Integer},
    ::Integer
)

end
