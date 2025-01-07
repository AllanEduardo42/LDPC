################################################################################
# Allan Eduardo Feitosa
# 14 Nov 2024
# Function to generate the exponent matrices for NR-LDPC encoding according to 
# the 3GPP 38.212 (Table 5.3.2.2)

include("BG.jl")

using LinearAlgebra

function
    make_exponent_matrix(n_BG,Zc,iLS)
    
    # 3GPP 38.212 Table 5.3.2.2
    if n_BG == 1
        n_rows_bg = 46; n_cols_bg = 68
        shift_coeffs_table = BG1()
    else
        n_rows_bg = 42; n_cols_bg = 52
        shift_coeffs_table = BG2()
    end
    
    
    #exponent matrix
    E_H  = -1 * ones(Int,n_rows_bg,n_cols_bg)
        
    # Loop all records in 
    n_rows_exp_mtx,_ = size(shift_coeffs_table)
    for i = 1:n_rows_exp_mtx 
        # row index and column index in exponent matrix
        row_idx = shift_coeffs_table[i,1]+1
        col_idx = shift_coeffs_table[i,2]+1
        
        # The corresponding shift V
        shift_coeff = shift_coeffs_table[i,iLS + 3]
        
        # Exponet matrix
        E_H[row_idx,col_idx] = rem(shift_coeff,Zc)
        
    end

    return E_H

end
    
    