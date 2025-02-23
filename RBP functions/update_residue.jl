################################################################################
# Allan Eduardo Feitosa
# 22 Feb 2024
# Function to find the maximum residue

# RBP
function 
    update_residue!(
        addressinv::Matrix{<:Integer},
        residues::Vector{<:AbstractFloat},
        l::Integer,
        x::AbstractFloat,
        Factors::Union{Matrix{<:AbstractFloat},Nothing}
    )

    @fastmath @inbounds x *= Factors[l]
    @inbounds residues[addressinv[l]] = x

end

# List-RBP
function 
    update_residue!(
        m::Integer,
        n::Integer,
        l::Integer,
        x::AbstractFloat,
        Factors::Union{Matrix{<:AbstractFloat},Nothing},
        listres1::Vector{<:AbstractFloat},
        listm1::Vector{<:Integer},
        listn1::Vector{<:Integer},
        listres2::Union{Vector{<:AbstractFloat},Nothing},
        listm2::Union{Vector{<:Integer},Nothing},
        listn2::Union{Vector{<:Integer},Nothing},
        listsize1::Integer,
        listsize2::Integer,
        new_listsize2::Integer,
        inlist::Matrix{<:Integer}
    )

    @fastmath @inbounds if inlist[l]  # if residue(m,n) is in the list
        # display("($m,$n) is on the list")
        inlist[l] = false   # remove from the list
        if listsize2 == 1
            new_listsize2 += 1
        end
        pos = 0
        for i = 1:listsize1
            if listm1[i] == m
                if listn1[i] == n
                    pos = i
                    break
                end
            end
        end
        if pos == 0
            throw(error("($m,$n) is on the list, but it's not registered."))
        end
        # remove from list
        remove_from_list!(l,listsize1,listres1,listm1,listn1,inlist,pos)
    end

    @fastmath @inbounds if x != 0.0
        x *= Factors[l]
        if listsize2 == 0
            update_list!(inlist,listres1,listm1,listn1,x,m,n,listsize1)
        elseif new_listsize2 == 1
            if x > listres2[1]
                listres2[1] = x
                listm2[1] = m
                listn2[1] = n
            end
        else
            update_list!(nothing,listres2,listm2,listn2,x,m,n,new_listsize2)
        end
    end

    return new_listsize2
    
end