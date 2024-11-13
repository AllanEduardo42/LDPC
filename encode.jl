################################################################################
# Allan Eduardo Feitosa
# 12 NOV 2024
# Function to generate parity bits

function 
    parity_bits(
        E_H::Matrix{<:Integer},
        u::Vector{Bool},
        Zc::Integer,
        bg::String)

    if bg == "1"
        I = 22
        J = 68
        M = 46
    else
        I = 10
        J = 52
        M = 42
    end

    E_A = E_H[1:4,1:I]
    E_C = E_H[5:end,1:I+4]

    a1 = circ_mult(E_A[1,:],I,Zc,u)
    a2 = circ_mult(E_A[2,:],I,Zc,u)
    a3 = circ_mult(E_A[3,:],I,Zc,u)
    a4 = circ_mult(E_A[4,:],I,Zc,u)

    p1 = a1 .⊻ a2 .⊻ a3 .⊻ a4

    if bg == "1"
        zp1 = circshift(p1,-1)
        p2 = circ_mult(E_H[1,:],I,Zc,u) .⊻ zp1
        p4 = circ_mult(E_H[4,:],I,Zc,u) .⊻ zp1
        p3 = circ_mult(E_H[3,:],I,Zc,u) .⊻ p4
        
    else
        p1 = circshift(p1,1)
        p2 = circ_mult(E_H[1,:],I,Zc,u) .⊻ p1
        p3 = circ_mult(E_H[2,:],I,Zc,u) .⊻ p2
        p4 = circ_mult(E_H[4,:],I,Zc,u) .⊻ p1
    end

    p = [p1;p2;p3;p4]

    c = [u;p] # K + 4*Zc

    q = zeros(Bool,(J-I-4)*Zc) # N - K - 4*Zc

    for m = 1:M-4

        q[1+Zc*(m-1):Zc*m] .= circ_mult(E_C[m,:],I+4,Zc,c)

    end

    return [p;q]

end

function circ_mult(e_A,I,Zc,x)

    y = zeros(Bool,Zc)

    for i = 1:I
        if e_A[i] != -1
            range = (1 + (i-1)*Zc) : (i*Zc)
            y .⊻= circshift(x[range],-e_A[i])
        end
    end

    return y
end