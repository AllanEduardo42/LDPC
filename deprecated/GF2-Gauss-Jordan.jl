## GF(2) Gauus-Jordan Elimination

function GF_GJ(H)

    M,N = size(H)

    A = copy(H)

    for i=1:M
        if A[i,i] != 1 && i < M
            p = i+1
            while p <= M && A[p,i] != 1
                p +=1
            end
            if p <= M
                row_i = A[i,:]
                row_p = A[p,:]
                A[i,:] = row_p
                A[p,:] = row_i
            end
        end
        for j=i+1:M
            if A[j,i] == 1
                A[j,:] = (A[j,:] + A[i,:]) .% 2
            end
        end
    end

    ## reduce

    B = copy(A)

    for i=M:-1:1
        p = i
        while p <= N && B[i,p] != 1
            p += 1
        end
        if p <= N
            for j=i-1:-1:1
                if B[j,p] == 1
                    B[j,:] = (B[j,:] + B[i,:]) .% 2
                end
            end
        end
    end
    
    return A, B

end


# test GF_GJ

M = 8
N = 12

H = rand([0,1],M,N)

A, B = GF_GJ(H)

