################################################################################
# Allan Eduardo Feitosa
# 19 set 2024
# Function to find the maximum residue for the RBP Sum-Product Algorithm

function 
    RBP_findmax_ΔLr!(
        _nodes::Vector{<:Integer},
        Lrn::Vector{<:AbstractFloat},
        Lq::AbstractVector{<:AbstractFloat},
        Lr::AbstractVector{<:AbstractFloat},
        check::Integer,
        nmax::Integer,
        max_residue::AbstractFloat,
        max_coords::Vector{<:Integer}
        )
    
    pLr = 1.0
    for node in _nodes
        @inbounds @fastmath Lrn[node] = tanh(0.5*Lq[node])
        @inbounds @fastmath pLr *= Lrn[node]
    end
    for node in _nodes
        if node != nmax
            @inbounds @fastmath x = pLr/Lrn[node]
            if abs(x) < 1
                @inbounds @fastmath ΔLr = 2*atanh(x) - Lr[node]
                if abs(ΔLr) > abs(max_residue)
                    max_residue = ΔLr                    
                    max_coords[1] = check
                    max_coords[2] = node
                end
            end
        end
    end

    return max_residue
end

function 
    calc_residues!(
        _nodes::Vector{<:Integer},
        Lrn::Vector{<:AbstractFloat},
        Lq::AbstractVector{<:AbstractFloat},
        Lr::AbstractVector{<:AbstractFloat},
        check::Integer,
        nmax::Integer,
        R::Matrix{<:AbstractFloat}
        )
    
    pLr = 1.0
    for node in _nodes
        @inbounds @fastmath Lrn[node] = tanh(0.5*Lq[node])
        @inbounds @fastmath pLr *= Lrn[node]
    end
    for node in _nodes
        if node != nmax
            @inbounds @fastmath x = pLr/Lrn[node]
            if abs(x) < 1
                @inbounds @fastmath R[check,node] = 2*atanh(x) - Lr[node]
            end
        end
    end

end