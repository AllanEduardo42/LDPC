include("BG.jl")
include("find_set_index_lift_size.jl")

using LinearAlgebra

function
    make_parity_check_matrix(n_BG,Zc)

    set_Idx = find_set_index_lift_size(Zc)

    return make_parity_check_matrix(n_BG,Zc,set_Idx)
    
end

function
    make_parity_check_matrix(n_BG,Zc,set_Idx)

    # [Ref] Nguyen, Tram Thi Bao, Tuy Nguyen Tan, and Hanho Lee. 
    #  "Efficient QC-LDPC encoder for 5G new radio." Electronics 8.6 (2019): 668.

    
    # 3GPP 38.212 Table 5.3.2.2(3)
    if n_BG == "1"
        n_rows_bg = 46; n_cols_bg = 68
        shift_coeffs_table = BG1()
    else
        n_rows_bg = 42; n_cols_bg = 52
        shift_coeffs_table = BG2()
    end
    
    
    #exponent matrix
    E_H  = -1 * ones(Int,n_rows_bg,n_cols_bg)
    H = zeros(Bool,n_rows_bg*Zc, n_cols_bg*Zc)
        
    # Loop all records in 
    n_rows_exp_mtx,_ = size(shift_coeffs_table)
    I_matrix = Matrix(I(Zc))
    for i = 1:n_rows_exp_mtx 
        # row index and column index in exponent matrix
        row_idx = shift_coeffs_table[i,1]+1
        col_idx = shift_coeffs_table[i,2]+1
        
        # The corresponding shift V
        shift_coeff = shift_coeffs_table[i,set_Idx + 3]
        
        # Exponet matrix
        E_H[row_idx,col_idx] = rem(shift_coeff,Zc)
        
        # check matrix
        row_range = Zc*(row_idx-1)+1 : Zc*row_idx
        col_range = Zc*(col_idx-1)+1 : Zc*col_idx
        H[row_range,col_range] = circshift(I_matrix,-rem(shift_coeff,Zc))
        
    end

    return BitMatrix(H), E_H

end
    
    