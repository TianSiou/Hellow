function [ mtx_deflate, flag ] = construct_deflated_mtx( lambda, ev_lambda, ev_lambda_inv, mtx_NME, sigma, LU_W, deflate_tol )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nu       = (lambda - 1)./(lambda + 1);
leng_nu  = length(nu);
idx_real = find(abs(imag(lambda)) <= 1.0e-10);
idx_pimg = find(abs(real(nu)) <= 1.0e-10);
idx_cmp  = setdiff(1:leng_nu, [idx_real; idx_pimg]);
m_1      = length(idx_real);
m_2      = length(idx_cmp); %leng_nu - length(idx_real);
m_3      = length(idx_pimg);
dim_m    = 2 * m_1 + 4 * m_2 + 2 * m_3;
Lambda_m = zeros(dim_m);
Lambda_d = zeros(2*m_1,1);
mtx_Xm   = zeros(size(ev_lambda,1), dim_m); 

% size(idx_pimg)
% size(idx_cmp)
% size(lambda)
% def_unit_ew = lambda([idx_pimg; idx_cmp'])
% figure(4)
% plot(def_unit_ew,'m^')

Lambda_d(1:m_1,1)         = real(nu(idx_real));
Lambda_d(m_1+1:2*m_1,1)   = -Lambda_d(1:m_1,1);
Lambda_m(1:2*m_1,1:2*m_1) = diag(Lambda_d);

for ii = 1:m_2
    Lambda_m(2*m_1+(ii-1)*2+1:2*m_1+2*ii,2*m_1+(ii-1)*2+1:2*m_1+2*ii) = ...
        [real(nu(idx_cmp(ii))) imag(nu(idx_cmp(ii))); -imag(nu(idx_cmp(ii))) real(nu(idx_cmp(ii)))]; 
end

kk = 2 * (m_1 + m_2);
Lambda_m(kk+1:kk+2*m_2,kk+1:kk+2*m_2) = - Lambda_m(2*m_1+1:kk,2*m_1+1:kk);

kk2 = 2 * m_1 + 4 * m_2;
for ii = 1:m_3
    Lambda_m(kk2+(ii-1)*2+1:kk2+2*ii,kk2+(ii-1)*2+1:kk2+2*ii) = ...
        [real(nu(idx_pimg(ii))) imag(nu(idx_pimg(ii))); -imag(nu(idx_pimg(ii))) real(nu(idx_pimg(ii)))]; 
end

%
% Construct othonormal basis Xm of the eigenspace
%
mtx_Xm(:,1:m_1)       = real(ev_lambda(:,idx_real));
mtx_Xm(:,m_1+1:2*m_1) = real(ev_lambda_inv(:,idx_real)); 

for ii = 1:m_2
    mtx_Xm(:,2*m_1+(ii-1)*2+1:2*m_1+2*ii) = ...
        [real(ev_lambda(:,idx_cmp(ii))) imag(ev_lambda(:,idx_cmp(ii)))]; 
end

for ii = 1:m_2
    mtx_Xm(:,kk+(ii-1)*2+1:kk+2*ii) = ...
        [real(ev_lambda_inv(:,idx_cmp(ii))) imag(ev_lambda_inv(:,idx_cmp(ii)))]; 
end

for ii = 1:m_3
    mtx_Xm(:,kk2+(ii-1)*2+1:kk2+2*ii) = ...
        [real(ev_lambda(:,idx_pimg(ii))) imag(ev_lambda(:,idx_pimg(ii)))]; 
end

% for ii = 1:length(lambda)
%    rsdl  = lambda(ii)^2 * (mtx_NME.mtx_A.' * ev_lambda(:,ii)) - lambda(ii) * (mtx_NME.mtx_Q * ev_lambda(:,ii)) + mtx_NME.mtx_A * ev_lambda(:,ii);
%    rsdl2 = lambda(ii)^(-2) * (mtx_NME.mtx_A.' * ev_lambda_inv(:,ii)) - lambda(ii)^(-1) * (mtx_NME.mtx_Q * ev_lambda_inv(:,ii)) + mtx_NME.mtx_A * ev_lambda_inv(:,ii);
%    fprintf('rsdl = %10.4e, rsdl_inv = %10.4e\n',norm(rsdl), norm(rsdl2));
% end
% 
 
% mtx_M1 = mtx_NME.mtx_A.' + mtx_NME.mtx_Q + mtx_NME.mtx_A;
% mtx_G = 2 * mtx_NME.mtx_A.' - 2 * mtx_NME.mtx_A;
% mtx_K = mtx_NME.mtx_A.' - mtx_NME.mtx_Q + mtx_NME.mtx_A;
%    
% Err = mtx_M1 * mtx_Xm * Lambda_m^2 + mtx_G * mtx_Xm * Lambda_m + mtx_K * mtx_Xm;
% % norm(Err(:,1:2*m_1),'fro')
% % norm(Err(:,1+2*m_1:2*m_1+2*m_2),'fro')
% fprintf('err of eigenmatrix = %11.4e \n',norm(Err,'fro'))
    
[mtx_Xm2, R] = qr(mtx_Xm,0);
% svd(R)
% svd(mtx_Xm(:,1+2*m_1:2*m_1+2*m_2))
% lambda(idx_cmp)

Lambda_m     = R * (Lambda_m / R);

tmp          = mtx_NME.mtx_A.' + mtx_NME.mtx_A;

% ATA_X   = tmp_mtx * mtx_Xm2;
% Q_X     = mtx_NME.mtx_Q * mtx_Xm2;
mtx_M1  = tmp + mtx_NME.mtx_Q;
mtx_G   = 2 * mtx_NME.mtx_A.' - 2 * mtx_NME.mtx_A;
mtx_K   = tmp - mtx_NME.mtx_Q;
    
Err  = mtx_M1 * mtx_Xm2 * Lambda_m^2 + mtx_G * mtx_Xm2 * Lambda_m + mtx_K * mtx_Xm2;
rsdl = norm(Err,'fro');
% 
% norm(Err(:,1:2*m_1),'fro')
% norm(Err(:,1+2*m_1:2*m_1+4*m_2),'fro')
fprintf('err of eigenmatrix = %11.4e \n',rsdl)

if ( rsdl < deflate_tol )
    flag                   = 0;
    mtx_deflate.mtx_Xm     = mtx_Xm2;
    mtx_deflate.Lambda_m   = Lambda_m;
    mtx_deflate.mtx_Vm     = (tmp + mtx_NME.mtx_Q) * mtx_deflate.mtx_Xm;
    Theta_inv              = mtx_deflate.mtx_Xm.' * mtx_deflate.mtx_Vm;
    Theta_inv              = (Theta_inv + Theta_inv.') / 2;
    mtx_deflate.mtx_Vm     = mtx_deflate.mtx_Vm * mtx_deflate.Lambda_m;
    mtx_deflate.mtx_Um     = (tmp - mtx_NME.mtx_Q) * mtx_deflate.mtx_Xm;
    mtx_deflate.mtx_VmTran = mtx_deflate.mtx_Vm.';
    mtx_deflate.mtx_UmTran = mtx_deflate.mtx_Um.';
    mtx_deflate.mtx_UTmVT  = mtx_deflate.mtx_UmTran - mtx_deflate.mtx_VmTran;
    mtx_deflate.mtx_UTpVT  = mtx_deflate.mtx_UmTran + mtx_deflate.mtx_VmTran;
    mtx_deflate.mtx_UpV    = mtx_deflate.mtx_UTpVT.';
    mtx_deflate.mtx_UmV    = mtx_deflate.mtx_UTmVT.';

    mtx_deflate.Phi_m      = Lambda_m.' * Theta_inv * Lambda_m;
    mtx_deflate.Phi_m      = (mtx_deflate.Phi_m + mtx_deflate.Phi_m.') / 2;

    %
    % Construct the matrices for Sherman-Morrison-Woodbury formula
    %
    tmp_mtx                = [ (sigma+1)^2*mtx_deflate.mtx_Um+(1-sigma^2)*mtx_deflate.mtx_Vm ...
                               (sigma^2-1)*mtx_deflate.mtx_Um-(sigma-1)^2*mtx_deflate.mtx_Vm ];
    mtx_deflate.s2UVTran   = [ (sigma+1)^2*mtx_deflate.mtx_Um.'+(1-sigma^2)*mtx_deflate.mtx_Vm.'; ...
                               (sigma^2-1)*mtx_deflate.mtx_Um.'-(sigma-1)^2*mtx_deflate.mtx_Vm.' ]; %tmp_mtx.';

    if ( abs(sigma) > 3.0e-3 )
        mtx_deflate.mtx_Wm = LU_W.Perm_amd * (LU_W.upper_U \ (LU_W.Low_L \ ...
            (LU_W.Perm_LU * tmp_mtx(LU_W.Perm_amd_vec,:))));

        mtx_deflate.mtx_Zm = LU_W.Perm_amd * (LU_W.Perm_LU_Tran * (LU_W.Low_L_Tran \ (LU_W.upper_U_Tran ...
            \ [mtx_deflate.mtx_Um(LU_W.Perm_amd_vec,:), mtx_deflate.mtx_Vm(LU_W.Perm_amd_vec,:)])));
    else
        mtx_deflate.mtx_Wm = LU_W.mtx_Qp \ tmp_mtx;
        mtx_deflate.mtx_Zm = LU_W.mtx_QpT \ [mtx_deflate.mtx_Um, mtx_deflate.mtx_Vm];
    end
    %
    % Construct small size matrix for Qp(lambda) in Sherman-Morrison-Woodbury formula
    %
    mtx_SMW                = 4 * [ mtx_deflate.Phi_m zeros(size(Lambda_m,1)); ...
        zeros(size(Lambda_m,1)) mtx_deflate.Phi_m ] + ...
        [ mtx_deflate.mtx_UmTran; mtx_deflate.mtx_VmTran ] * mtx_deflate.mtx_Wm;

    [mtx_deflate.LU_L, mtx_deflate.LU_U, mtx_deflate.LU_P] = lu(mtx_SMW);
    %
    % Construct small size matrix for Qp(lambda)^T in Sherman-Morrison-Woodbury formula
    %
    mtx_SMWTran           = 4 * [ mtx_deflate.Phi_m zeros(size(Lambda_m,1)); ...
        zeros(size(Lambda_m,1)) mtx_deflate.Phi_m ] + mtx_deflate.s2UVTran * mtx_deflate.mtx_Zm;

    [mtx_deflate.TranLU_L, mtx_deflate.TranLU_U, mtx_deflate.TranLU_P] = lu(mtx_SMWTran);
else
    flag        = 1;
    mtx_deflate = [];
end
% tmp_mtx = mtx_NME.mtx_A.' + mtx_NME.mtx_A;
% % ATA_X   = tmp_mtx * mtx_Xm2;
% % Q_X     = mtx_NME.mtx_Q * mtx_Xm2;
% mtx_M1  = tmp_mtx + mtx_NME.mtx_Q;
% mtx_G   = 2 * mtx_NME.mtx_A.' - 2 * mtx_NME.mtx_A;
% mtx_K   = tmp_mtx - mtx_NME.mtx_Q;
%     
% Err = mtx_M1 * mtx_Xm2 * Lambda_m^2 + mtx_G * mtx_Xm2 * Lambda_m + mtx_K * mtx_Xm2;
% % size(Err)
% % 2*m_1+4*m_2+2*m_3
% norm(Err(:,1:2*m_1),'fro')
% norm(Err(:,1+2*m_1:2*m_1+4*m_2),'fro')
% fprintf('err of eigenmatrix = %11.4e \n',norm(Err,'fro'))

end

