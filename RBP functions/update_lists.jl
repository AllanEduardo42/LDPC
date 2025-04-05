################################################################################
# Allan Eduardo Feitosa
# 2 Mar 2025
# Update the list of residues for the List-RBP algorithm

include("add_to_list.jl")

function
    update_list2!(
        residues::Vector{<:AbstractFloat},
        coords::Matrix{<:Integer},
        localresidues::Vector{<:AbstractFloat},
        localcoords::Matrix{<:Integer},
        listsizes::Vector{<:Integer},
        rbpmatrix::Matrix{<:Integer},
        li::Integer,
        m::Integer,
        n::Integer,
        residue::AbstractFloat    
    )
    
    @inbounds begin
        if rbpmatrix[li]  # if residue(m,n) is in the list
            # display("($m,$n) is on the list")
            pos = 0
            for i = 1:listsizes[1]
                if coords[3,i] == li
                    pos = i
                    break
                end
            end
            if pos == 0
                throw(error("($(coords[1,i]),$(coords[1,i])) is registered as being on the list, but it's not."))
            end
            # remove from list
            rbpmatrix[li] = false
            remove_maxresidue!(li,listsizes[1],residues,coords,rbpmatrix,pos)
        end                
        if residue != 0.0
            add_to_list!(nothing,localresidues,localcoords,residue,li,m,n,listsizes[2])
        end
    end

end

function
    update_list2!(
        residues::Vector{<:AbstractFloat},
        ::Matrix{<:Integer},
        ::Nothing,
        ::Nothing,
        ::Vector{<:Integer},
        rbpmatrix::Matrix{<:Integer},
        li::Integer,
        m::Integer,
        n::Integer,
        residue::AbstractFloat    
    )

    @inbounds residues[rbpmatrix[li]] = residue

    return 0
end

function update_list1!(
    residues::Vector{<:AbstractFloat},
    coords::Matrix{<:Integer},
    localresidues::Vector{<:AbstractFloat},
    localcoords::Matrix{<:Integer},
    listsizes::Vector{<:Integer},
    rbpmatrix::Matrix{<:Integer}
)
    
    @inbounds begin
        for i = listsizes[3]:-1:listsizes[4]
            li = coords[3,i]
            if li ≠ 0
                rbpmatrix[li] = false
                residues[i] = 0.0
            end
        end
        for i in 1:listsizes[2]
            m = localcoords[1,i]
            n = localcoords[2,i]
            li = localcoords[3,i]
            add_to_list!(rbpmatrix,residues,coords,localresidues[i],li,m,n,listsizes[1])
        end
        # clear list 2
        localresidues .*= 0.0
        localcoords .*= 0
    end
end

function update_list1!(
    ::Vector{<:AbstractFloat},
    ::Matrix{<:Integer},
    ::Nothing,
    ::Nothing,
    ::Vector{<:Integer},
    ::Matrix{<:Integer}
)

end