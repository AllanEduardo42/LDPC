################################################################################
# Allan Eduardo Feitosa
# 2 Mar 2025
# List-RBP auxiliary functions

# Add the residual to the list
function add_to_list!(
    inlist::Union{Matrix{Bool},Nothing},
    list::Vector{Float64},
    coords::Matrix{Int},
    residual::Float64,
    li::Int,
    ci::Int,
    vj::Int,
    listsize::Int
)

    @fastmath @inbounds if residual > list[listsize]

        if residual ≥ list[1]
            i = 1
        else
            d = listsize >> 1
            i = d
            while d > 1
                d >>= 1
                if residual ≥ list[i]
                    i -= d
                else
                    i += d
                end
            end
            if residual < list[i]
                i += 1
            end
        end

        update_inlist!(inlist,coords,li,listsize)

        for j=listsize:-1:i+1
            list[j] = list[j-1]
            coords[1,j] = coords[1,j-1]
            coords[2,j] = coords[2,j-1]
            coords[3,j] = coords[3,j-1]
        end


        coords[1,i] = ci
        coords[2,i] = vj
        coords[3,i] = li
        list[i] = residual        
    end
end

# update inlist matrix
function update_inlist!(
    inlist::Matrix{Bool},
    coords::Matrix{Int},
    li::Int,
    listsize::Int
)

     begin
        last = coords[3,listsize]
        if last ≠ 0
            inlist[last] = false
        end
        inlist[li] = true
    end
end

function update_inlist!(
    ::Nothing,
    ::Matrix{Int},
    ::Int,
    ::Int
)

end

# Find position (index) in the list
function find_list_pos(
    li::Int,
    listsize::Int,
    coords::Matrix{<:Int},
    ci::Int,
    vj::Int
)

    pos = 0
    @inbounds for i = 1:listsize
        if coords[3,i] == li
            pos = i
            break
        end
    end
    if pos == 0
        throw(error("($ci,$vj) is registered as being on the list, but it's not."))
    end

    return pos

end

# Initializes main list
function init_list!(
    V2C::Matrix{Float64},
    Nc::Vector{Vector{Int}},
    phi::Union{Vector{Float64},Nothing},
    msum_factor::Union{Float64,Nothing},
    newC2V::Matrix{Float64},   
    inlist::Matrix{Bool},
    Residuals::Matrix{Float64},
    list::Vector{Float64},
    coords::Matrix{Int},
    listsize::Int;
    listsize2=listsize
)

    @fastmath @inbounds for ci in eachindex(Nc)
        Nci = Nc[ci]
        for vj in Nci
            li = LinearIndices(newC2V)[ci,vj]
            newlr = calc_C2V(Nci,ci,vj,V2C,msum_factor)
            newC2V[li] = newlr
            residual = abs(newlr)
            Residuals[li] = residual
            add_to_list!(inlist,list,coords,residual,li,ci,vj,listsize)
        end
    end
end

# removes item from the list
function remove_from_list!(
    li::Int,
    listsize::Int,
    list::Vector{Float64},
    coords::Matrix{Int},
    inlist::Matrix{Bool},
    pos::Int
)

    @inbounds inlist[li] = false

    # update the list
    @inbounds for i in pos:listsize
        list[i] = list[i+1]
        coords[1,i] = coords[1,i+1]
        coords[2,i] = coords[2,i+1]
        coords[3,i] = coords[3,i+1]
    end    
end

# Update main list
function update_main_list!(
    list::Vector{Float64},
    coords::Matrix{Int},
    local_list::Vector{Float64},
    local_coords::Matrix{Int},
    listsize::Int,
    listsize2::Int,
    inlist::Matrix{Bool}
)
    
    @inbounds begin
        pos = listsize-1
        for i=listsize2:-1:1
            if local_coords[3,i] != 0
                pos = listsize - i + 1
                break
            end
        end
        for i = pos:listsize-1
            li = coords[3,i]
            if li ≠ 0
                inlist[li] = false
                list[i] = 0.0
            else
                break
            end
        end
        for i in 1:listsize2
            m = local_coords[1,i]
            n = local_coords[2,i]
            li = local_coords[3,i]
            add_to_list!(inlist,list,coords,local_list[i],li,m,n,listsize)
        end
        # clear list 2
        local_list .*= 0.0
        local_coords .*= 0
    end
end


