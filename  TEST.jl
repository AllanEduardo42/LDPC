# TEST

function residue1(A,B,C)

    for m in axes(B,1)
        for n in axes(B,2)
            @fastmath @inbounds A[m,n] = abs(B[m,n] - C[m,n])
        end
    end
end

function residue2(A,B,C)
    
    for n in axes(B,2)
        colB = B[:,n]
        colC = C[:,n]
        for m in axes(B,1)
            @fastmath @inbounds A[m,n] = abs(colB[m] - colC[m])
        end
    end
end

function residue3(A,B,C)
    
    for m in axes(B,1)
        rowB = B[m,:]
        rowC = C[m,:]
        for n in axes(B,2)
            @fastmath @inbounds A[m,n] = abs(rowB[n] - rowC[n])
        end
    end
end

function residue4(A,B,C)
    
    for n in axes(B,2)
        colB = view(B,:,n)
        colC = view(C,:,n)
        for m in axes(B,1)
            @fastmath @inbounds A[m,n] = abs(colB[m] - colC[m])
        end
    end
end

function residue5(A,B,C)

    for n in axes(B,2)
        for m in axes(B,1)
            @fastmath @inbounds A[m,n] = abs(B[m,n] - C[m,n])
        end
    end
end


M = 1024
N = 2048
A = zeros(M,N)
B = randn(M,N)
C = randn(M,N)

@time residue5(A,B,C)
@time residue1(A,B,C)
@time residue4(A,B,C)
@time residue2(A,B,C)
@time residue3(A,B,C)


@time residue5(A,B,C)
@time residue1(A,B,C)
@time residue4(A,B,C)
@time residue2(A,B,C)
@time residue3(A,B,C)


