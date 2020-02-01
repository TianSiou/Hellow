clear all
clc

D11 = rand(10,1); D12 = rand(10,1); D13 = rand(10,1);
D21 = rand(10,1); D22 = rand(10,1); D23 = rand(10,1);
D31 = rand(10,1); D32 = rand(10,1); D33 = rand(10,1);

[D, invD] = Matrix_Inverse_BlockMatrix(D11,D12,D13,D21,D22,D23,D31,D32,D33);