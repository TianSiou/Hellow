function [ sol ] = Defl_mtx_vec_A( vec, transp_flag, mtx_NME, mtx_deflate )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if strcmp(transp_flag,'transp')
    sol = 4 * (mtx_NME.mtx_A.' * vec) + mtx_deflate.mtx_UmV * (mtx_deflate.Phi_m \ (mtx_deflate.mtx_UTpVT * vec));
else
    sol = 4 * (mtx_NME.mtx_A * vec) + mtx_deflate.mtx_UpV * (mtx_deflate.Phi_m \ (mtx_deflate.mtx_UTmVT * vec));
end

end

