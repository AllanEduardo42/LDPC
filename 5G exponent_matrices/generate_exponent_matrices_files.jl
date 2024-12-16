################################################################################
# Allan Eduardo Feitosa
# 14 Nov 2024
# Routine to generate the files of the exponent matrices of NR-LDPC encoding,
# according to the 3GPP 38.212 (Tables 5.3.2.1 and 5.3.2.2)

include("liftSizeMtx.jl")
include("make_exponent_matrix_matrix.jl")

using DelimitedFiles

for i in axes(liftSizeMtx,1)
    iLS = i - 1
    for j in axes(liftSizeMtx,2)
        Zc = liftSizeMtx[i,j]
        if Zc ≠ 0
            for bg = 1:2
                E_H = make_exponent_matrix(bg,Zc,iLS)
                open("./exponent_matrices/EM_$(bg)_$(iLS)_$(Zc).txt","w") do io
                    writedlm(io,E_H)
                end
            end
        end
    end
end