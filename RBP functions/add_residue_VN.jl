################################################################################
# Allan Eduardo Feitosa
# 
#

# List-VN-RBP
function add_residue_VN!(
        inlist::Union{Vector{Bool},Nothing},
        alpha::Vector{<:AbstractFloat},
        coords::Vector{<:Integer},
        residue::AbstractFloat,
        vj::Integer,
        listsize::Integer
    )

    @fastmath @inbounds if residue > alpha[listsize]

        if listsize > 1
            if residue ≥ alpha[1]
                i = 1
            else
                d = listsize >> 1
                i = d
                while d > 1
                    d >>= 1
                    if residue ≥ alpha[i]
                        i -= d
                    else
                        i += d
                    end
                end
                if residue < alpha[i]
                    i += 1
                end
            end

            update_inlist_VN!(inlist,coords,vj)

            for j=listsize:-1:i+1
                alpha[j] = alpha[j-1]
                coords[j] = coords[j-1]
            end
        else
            i = 1
        end

        coords[i] = vj
        alpha[i] = residue        
    end
end

# auxiliary function
function update_inlist_VN!(
    inlist::Vector{Bool},
    coords::Vector{<:Integer},
    vj::Integer
)

    @inbounds begin
        last = coords[end-1]
        if last ≠ 0
            inlist[last] = false
        end
        inlist[vj] = true
    end
end

function update_inlist_VN!(
    ::Nothing,
    ::Vector{<:Integer},
    ::Integer
)

end
