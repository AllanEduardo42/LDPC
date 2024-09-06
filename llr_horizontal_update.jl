################################################################################
# Allan Eduardo Feitosa
# 27 ago 2024
# Horizontal update of the LLR based Sum-Product Algorithm

using SparseArrays

# TANH
function llr_horizontal_update!(Lr::Matrix{Float64},
                                Lq::Matrix{Float64},
                                indices_n::Vector{Vector{Int64}},
                                x::Nothing,
                                y::Nothing,
                                z::Nothing)

    m = 0
    @inbounds for indices in indices_n
        m += 1
        @inbounds for n in indices
            pLr = 1.0
            @inbounds for nn in indices
                if nn != n
                    @fastmath pLr *= tanh(0.5*Lq[nn,m])
                end
            end
            @fastmath Lr[m,n] = 2*atanh(pLr)
        end
    end
end

# APPROX
function llr_horizontal_update!(Lr::Matrix{Float64},                           
                                Lq::Matrix{Float64},
                                indices_n::Vector{Vector{Int64}},
                                sn::Vector{Int64},
                                x::Nothing,
                                y::Nothing)    

    m = 0
    @inbounds for indices in indices_n
        m += 1
        s = 1
        minL = 0.0
        minL2 = 0.0
        idx = 0
        @inbounds for n in indices
            Lrn = Lq[n,m]
            if Lrn > 0
                sn[n] = 1
            else
                sn[n] = -1
                s *= sn[n]
                Lrn *= -1
            end
            if minL == 0.0
                idx = n
                minL = Lrn
            elseif minL2 == 0
                if Lrn < minL
                    idx = n
                    minL2 = minL
                    minL = Lrn
                else
                    minL2 = Lrn
                end
            else
                if Lrn < minL
                    idx = n
                    minL2 = minL
                    minL = Lrn
                elseif Lrn < minL2
                    minL2 = Lrn
                end
            end
        end
        @inbounds for n in indices
            if n == idx
                if sn[n] != s
                    Lr[m,n] = -minL2
                else
                    Lr[m,n] = minL2
                end
            else
                if sn[n] != s
                    Lr[m,n] = -minL
                else
                    Lr[m,n] = minL
                end
            end
        end  
    end
end

# ALT
function llr_horizontal_update!(Lr::Matrix{Float64},                            
                                Lq::Matrix{Float64},
                                indices_n::Vector{Vector{Int64}},
                                sn::Vector{Int64},
                                Lrn::Vector{Float64},
                                x::Nothing)

    m = 0
    @inbounds for indices in indices_n
        m += 1
        sum = 0.0
        s = 1 
        @inbounds for n in indices
            ####### what it does: #######
            # Lrn[n] = ϕ(abs(Lq[n,m]))
            # sum += Lrn[n]
            # sn[n] = sign(Lq[n,m])
            # s *= sn[n]
            #############################
            if Lq[n,m] >= 0
                Lrn[n] = ϕ(Lq[n,m])
                sn[n] = 1
            else
                Lrn[n] = ϕ(-Lq[n,m])
                sn[n] = -1
                s *= sn[n]
            end
            @fastmath sum += Lrn[n]
        end
        @inbounds for n in indices
            ############### what it does ###############
            # Lr[m,n] = (s*sn[n])*(ϕ(abs(sum - Lrn[n])))
            ############################################
            if sn[n] != s
                @fastmath Lr[m,n] = -ϕ(abs(sum - Lrn[n]))
            else
                @fastmath Lr[m,n] = ϕ(abs(sum - Lrn[n]))
            end
        end    
    end
end

# TABLE
function llr_horizontal_update!(Lr::Matrix{Float64},                            
                                Lq::Matrix{Float64},
                                indices_n::Vector{Vector{Int64}},
                                sn::Vector{Int64},
                                Lrn::Vector{Float64},
                                phi::Vector{Float64})

    m = 0
    @inbounds for indices in indices_n
        m += 1
        sum = 0.0
        s = 1 
        @inbounds for n in indices
            if Lq[n,m] >= 0
                Lrn[n] = phi[get_index(Lq[n,m])]
                sn[n] = 1
            else
                Lrn[n] = phi[get_index(-Lq[n,m])]
                sn[n] = -1
                s *= sn[n]
            end
            @fastmath sum += Lrn[n]
        end
        @inbounds for n in indices
            if sn[n] != s
                @fastmath Lr[m,n] = -phi[get_index(abs(sum - Lrn[n]))]
            else
                @fastmath Lr[m,n] = phi[get_index(abs(sum - Lrn[n]))]
            end
        end    
    end
end

# TANH
function llr_horizontal_update!(Lr::SparseMatrixCSC{Float64, Int64},
                                Lq::SparseMatrixCSC{Float64, Int64},
                                indices_n::Vector{Vector{Int64}},
                                x::Nothing,
                                y::Nothing,
                                z::Nothing)

    m = 0
    @inbounds for indices in indices_n
        m += 1
        @inbounds for n in indices
            pLr = 1.0
            @inbounds for nn in indices
                if nn != n
                    @fastmath pLr *= tanh(0.5*Lq[nn,m])
                end
            end
            @fastmath Lr[m,n] = 2*atanh(pLr)
        end
    end
end

# APPROX
function llr_horizontal_update!(Lr::SparseMatrixCSC{Float64, Int64},                           
                                Lq::SparseMatrixCSC{Float64, Int64},
                                indices_n::Vector{Vector{Int64}},
                                sn::Vector{Int64},
                                x::Nothing,
                                y::Nothing)    

    m = 0
    @inbounds for indices in indices_n
        m += 1
        s = 1
        minL = 0.0
        minL2 = 0.0
        idx = 0
        @inbounds for n in indices
            Lrn = Lq[n,m]
            if Lrn > 0
                sn[n] = 1
            else
                sn[n] = -1
                s *= sn[n]
                Lrn *= -1
            end
            if minL == 0.0
                idx = n
                minL = Lrn
            elseif minL2 == 0
                if Lrn < minL
                    idx = n
                    minL2 = minL
                    minL = Lrn
                else
                    minL2 = Lrn
                end
            else
                if Lrn < minL
                    idx = n
                    minL2 = minL
                    minL = Lrn
                elseif Lrn < minL2
                    minL2 = Lrn
                end
            end
        end
        @inbounds for n in indices
            if n == idx
                if sn[n] != s
                    Lr[m,n] = -minL2
                else
                    Lr[m,n] = minL2
                end
            else
                if sn[n] != s
                    Lr[m,n] = -minL
                else
                    Lr[m,n] = minL
                end
            end
        end  
    end
end

# ALT
function llr_horizontal_update!(Lr::SparseMatrixCSC{Float64, Int64},                            
                                Lq::SparseMatrixCSC{Float64, Int64},
                                indices_n::Vector{Vector{Int64}},
                                sn::Vector{Int64},
                                Lrn::Vector{Float64},
                                x::Nothing)

    m = 0
    @inbounds for indices in indices_n
        m += 1
        sum = 0.0
        s = 1 
        @inbounds for n in indices
            ####### what it does: #######
            # Lrn[n] = ϕ(abs(Lq[n,m]))
            # sum += Lrn[n]
            # sn[n] = sign(Lq[n,m])
            # s *= sn[n]
            #############################
            if Lq[n,m] >= 0
                Lrn[n] = ϕ(Lq[n,m])
                sn[n] = 1
            else
                Lrn[n] = ϕ(-Lq[n,m])
                sn[n] = -1
                s *= sn[n]
            end
            @fastmath sum += Lrn[n]
        end
        @inbounds for n in indices
            ############### what it does ###############
            # Lr[m,n] = (s*sn[n])*(ϕ(abs(sum - Lrn[n])))
            ############################################
            if sn[n] != s
                @fastmath Lr[m,n] = -ϕ(abs(sum - Lrn[n]))
            else
                @fastmath Lr[m,n] = ϕ(abs(sum - Lrn[n]))
            end
        end    
    end
end

# TABLE
function llr_horizontal_update!(Lr::SparseMatrixCSC{Float64, Int64},                            
                                Lq::SparseMatrixCSC{Float64, Int64},
                                indices_n::Vector{Vector{Int64}},
                                sn::Vector{Int64},
                                Lrn::Vector{Float64},
                                phi::Vector{Float64})

    m = 0
    @inbounds for indices in indices_n
        m += 1
        sum = 0.0
        s = 1 
        @inbounds for n in indices
            if Lq[n,m] >= 0
                Lrn[n] = phi[get_index(Lq[n,m])]
                sn[n] = 1
            else
                Lrn[n] = phi[get_index(-Lq[n,m])]
                sn[n] = -1
                s *= sn[n]
            end
            @fastmath sum += Lrn[n]
        end
        @inbounds for n in indices
            if sn[n] != s
                @fastmath Lr[m,n] = -phi[get_index(abs(sum - Lrn[n]))]
            else
                @fastmath Lr[m,n] = phi[get_index(abs(sum - Lrn[n]))]
            end
        end    
    end
end