function add_to_list!(
        rbpmatrix::Union{Matrix{Bool},Nothing},
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

            update_inlist!(rbpmatrix,coords,li)

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

function add_to_list!(
    rbpmatrix::Matrix{<:Integer},
    residues::Vector{<:AbstractFloat},
    ::Matrix{<:Integer},
    residue::AbstractFloat,
    li::Integer,
    ::Integer,
    ::Integer,
    ::Integer
)

    @inbounds residues[rbpmatrix[li]] = residue

end

function update_inlist!(
    rbpmatrix::Matrix{Bool},
    coords::Matrix{<:Integer},
    li::Integer
)

    @inbounds begin
        last = coords[3,end-1]
        if last ≠ 0
            rbpmatrix[last] = false
        end
        rbpmatrix[li] = true

    end

end

function update_inlist!(
    ::Nothing,
    ::Matrix{<:Integer},
    ::Integer
)

end
