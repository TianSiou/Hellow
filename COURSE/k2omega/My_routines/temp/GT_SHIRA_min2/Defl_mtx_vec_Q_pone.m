function [ sol ] = Defl_mtx_vec_Q_pone( vec, mtx_NME, mtx_deflate )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

sol = 4 * (mtx_NME.mtx_Q * vec) + 2 * ([mtx_deflate.mtx_Um, mtx_deflate.mtx_Vm] * ...
    ([ mtx_deflate.Phi_m\(mtx_deflate.mtx_UmTran * vec); mtx_deflate.Phi_m\(mtx_deflate.mtx_VmTran * vec) ]));

end

