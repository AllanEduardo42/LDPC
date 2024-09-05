### Moreira's example using Mckay method:

H = [0 1 0 1 0 1 1 1 0 0 0 1;
     1 0 1 1 0 0 0 0 1 0 0 0;
     0 1 0 0 1 0 1 0 0 0 0 1;
     1 0 0 1 0 0 0 0 0 1 1 0;
     0 0 1 0 1 1 0 0 0 1 0 0;
     1 0 1 0 0 0 1 1 0 0 1 0;
     0 1 0 0 0 1 0 1 1 1 0 0;
     0 0 0 0 1 0 0 0 1 0 1 1]

MM,NN = size(H)
K = NN - MM

A = H[:,1:MM]
B = H[:,MM+1:NN]

P = abs.(Int64.(A\B .% 2))

G = [P; I(K)]

girth = 4