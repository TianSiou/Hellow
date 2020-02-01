function [ lambda_tmp, residual, R, H, Y, Z, nrm_Err_N ] = estimate_by_GTSHIRA( mtx_A, mtx_Q, sigma, mtx_NME, solve_LS_PQEP, n, nmax )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% W                  = mtx_NME.mtx_A - sigma*mtx_NME.mtx_Q + (sigma^2)*(mtx_NME.mtx_A.');
% LU_W.Perm_amd_vec  = amd(W);
% W_reorder          = W(LU_W.Perm_amd_vec, LU_W.Perm_amd_vec);
% LU_W.Perm_amd      = sparse(LU_W.Perm_amd_vec,1:n,ones(n,1));
%         
% [LU_W.Low_L,LU_W.upper_U,LU_W.Perm_LU] = lu(W_reorder);
%         
% LU_W.Perm_LU_Tran  = LU_W.Perm_LU.';
% LU_W.Low_L_Tran    = LU_W.Low_L.';
% LU_W.upper_U_Tran  = LU_W.upper_U.'; 
% LU_W.sigma         = sigma;
%         
% if ( abs(sigma) <= 3.0e-3 )
%     LU_W.mtx_Qp  = W;
%     LU_W.mtx_QpT = W.';
% end
        
%fprintf('GTSHIRA with sigma = %11.4e + (%11.4e)*1i \n', real(sigma), imag(sigma));
  
s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);
xini   = randn(2*n,1) + 1i * randn(2*n,1);
xini   = xini / norm(xini);
    
n2     = 2 * n; 
%nmax   = 80;
R      = zeros(nmax);
H      = zeros(nmax+1,nmax);
Y      = zeros(n2,nmax+1);
Z      = zeros(n2,nmax); 

Y(:,1) = xini; 
Y(:,1) = Y(:,1) / norm(Y(:,1));

loopstart  = 1; loopend = nmax; 
[Y(:,1:loopend+1), Z(:,1:loopend), H(1:loopend+1,1:loopend), R(1:loopend,1:loopend)] =  ... 
     GTArnoldi_step(n, mtx_NME.mtx_A, mtx_NME.mtx_Q, solve_LS_PQEP, Y(:,1:loopstart), Z(:,1:loopstart-1), H(1:loopstart,1:loopstart-1), ...
     R(1:loopstart-1,1:loopstart-1), sigma, loopstart, loopend);

%fprintf('\n Error check for Krylov subspace in estimate_by_GTSHIRA \n'); 
% KZ    = mtx_K_prod_vec( mtx_A, Z(:,1:loopend), sigma );
% Err_K = KZ - Y(:,1:loopend) * H(1:loopend,1:loopend) - ...
%         H(1+loopend,loopend) * Y(:,loopend+1) * [ zeros(1,loopend-1) 1];
% nrm_Err_K = norm(Err_K,inf);  
% fprintf('Error of Arnoldi decomp with K = %11.4e\n',nrm_Err_K);

NZ        = mtx_N_prod_vec( mtx_A, mtx_Q, Z(:,1:loopend), sigma );
Err_N     = NZ - Y(:,1:loopend) * R(1:loopend,1:loopend);
nrm_Err_N = norm(Err_N,inf);
fprintf('Error of Arnoldi decomp with N = %11.4e in estimate_by_GTSHIRA \n\n',nrm_Err_N);


[REV,REW] = eig(H(1:nmax,1:nmax), R(1:nmax,1:nmax));
residual  = abs(H(nmax+1,nmax)) * abs(REV(nmax,:)); 
ew        = diag(REW); 

if (  nrm_Err_N > 1.0e-8 )
    if ( nrm_Err_N < 1 )
        tolerance =  1.0e-3 * abs(sigma);
    else
        tolerance =  1.0e-3 * abs(sigma) / nrm_Err_N;
    end
elseif ( abs(sigma) > 1.0e-3 )
    tolerance =  1.0e-2 * abs(sigma)^(0.3); %abs(sigma)^(0.3);
elseif ( mod(angle(sigma), 2*pi) < 165/180*pi )
    tolerance =  1.0e-2 * abs(sigma)^(0.5);
else
    tolerance =  1.0e-2 * abs(sigma)^(0.8);
end

residual  = residual./abs(ew).';
idx_conv  = find(residual <= tolerance);
if ( length(idx_conv) < 10 )
    idx_conv  = find(residual <= 5.0e-2);
end
 
if ( ~isempty(idx_conv) )
    ew          = ew(idx_conv);
    residual    = residual(idx_conv); 
end

%
% Compute the eigenvalues lambda_vec and lambda_vec.^(-1) of the symplectic pencil (M, L)
%
lambda_vec = sigma + (1.0/sigma) + 1./ew;
tmp        = lambda_vec.^2 - 4;
lambda_vec = (lambda_vec + tmp.^(0.5))/2;
    
lambda_tmp = lambda_vec;
no_ew      = length(lambda_vec);
tol        = 1.0e-7;
for jj = 1:no_ew 
    lambda = lambda_vec(jj,1);
    if ( abs(lambda) > 1+tol )
        lambda_tmp(jj, 1) = 1/lambda; 
    elseif ( abs(lambda) < 1-tol )
        lambda_tmp(jj, 1) = lambda  ; 
    else
        if ( imag(lambda) > 0 )
            lambda_tmp(jj, 1) = lambda  ; 
        else
            lambda_tmp(jj, 1) = 1/lambda; 
        end
    end 
end
    
[~,idx_sort] = sort(real(lambda_tmp));
lambda_tmp   = lambda_tmp(idx_sort);
residual     = residual(idx_sort);
        
tmp          = abs(lambda_tmp + 1./lambda_tmp - sigma - 1 / sigma);
[~,idx]      = sort(tmp);
lambda_tmp   = lambda_tmp(idx(1:min(60,no_ew)));
residual     = residual(idx(1:min(60,no_ew)));
 
nstp                    = length(lambda_tmp);
tmp_H                   = H(loopend+1,loopend);
vec_y                   = Y(:,loopend+1);
                      
[new_H, new_R, Y_new, Z_new, alpha, new_nstp] = Schur_restart(nmax, nstp, H, R, Y, Z);
    
R = new_R;
H = [ new_H; zeros(1,new_nstp-1) alpha*tmp_H ]; 
Y = [Y_new  vec_y]; 
Z = Z_new;
% R(1:new_nstp,1:new_nstp) = new_R;
% H(1:new_nstp,1:new_nstp) = new_H;
% H(1+new_nstp,1:new_nstp) = [ zeros(1,new_nstp-1) alpha*tmp_H ];
% Y(:,1:new_nstp)          = Y_new;
% Y(:,1+new_nstp)          = vec_y;
% Z(:,1:new_nstp)          = Z_new;
nstp                     = new_nstp;
  
if (  nrm_Err_N > 1.0e-8 )
    NZ        = mtx_N_prod_vec( mtx_A, mtx_Q, Z, sigma );
    Err_N     = NZ - Y(:,1:nstp) * R;
    nrm_Err_N = norm(Err_N,inf);
    fprintf('Restarting, Error of Arnoldi decomp with N = %11.4e\n\n',nrm_Err_N);
end

end

