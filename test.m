clear all
clc
addpath('C:\Users\User\Desktop\FAME_m_g071708 - ½Æ»s (2)\Sira');
N = 1000;
 stop_tolerance = 1.0e-8; %1.0e8 * eps / scale; 
        target_type    = 'RGTR'; %'RGTC'; 

        no_restart     = 35; %35;
%         mtxdim         = 16 * n;
        CorEq          = 'SIRA'; %'JD'; %'SIRA';
        LSinfo.solver  = 'minres';%'bicgstabl'; %'bicgstabl'; %'minres';
        LSinfo.precond = 'no'; %'yes';
        flag_LU = 0;
        initial_V   = randn(3*N,1); % ones(dim,1);
        initial_V   = initial_V / norm(initial_V);
        fid = 0;
% val = ones(3*N,1);
val = [];
i_idx = []; j_idx = []; c = 1;
for i = 1:3*N
    val   = [val,c];
    i_idx = [i_idx,c];
    j_idx = [j_idx,c];
    c = c + 1;
end
% val = 1 ./ val;
A = sparse(i_idx, j_idx, val, 3*N, 3*N);
[ew, ev] = GEP_AB_Herm_JDSIRA_Driver (A,...
                         @(x) x, 3*N, no_restart, ...
            15, stop_tolerance, initial_V, 49.5, fid, target_type, CorEq, LSinfo); % @(x)solve_Minv_b( x, LU_precond ));     