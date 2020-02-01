function [new_H, new_R, Y_new, Z_new, alpha, new_nstp] = Schur_restart(nmax, nstp, H, R, Y, Z)

% global mtx_G mtx_M2 mtx_F mtx_M1 shift

[RH, RR, Qs, Zs] = qz(H(1:nmax,1:nmax), R(1:nmax,1:nmax));
alpha            = diag(RH);
beta             = diag(RR);
lambda           = alpha./beta;

clus_order       = zeros(nmax,1);
weight           = zeros(nmax,1);

[ew, sort_idx]   = sort(abs(lambda),'descend');
s                = nmax;
weight(1)        = s;

%
% How many |lambda| in ew are equal to ew(1), said 2. 
% Save 2 to multiplicity(1)
%
kk                = 1;
multiplicity      = zeros(nmax,1);
multiplicity(1,1) = 1;

for i = 2:nmax
    if ew(i) == ew(i-1)
        weight(i)          = s;
        multiplicity(kk,1) = multiplicity(kk,1) + 1;
    else
        s                  = s-1;
        weight(i)          = s;
        kk                 = kk + 1;
        multiplicity(kk,1) = 1;
    end
end
clus_order(sort_idx) = weight;
ct                   = cumsum(multiplicity(1:kk,1));
new_size             = find(ct >= nstp,1);
new_nstp             = ct(new_size);

[newRH, newRR, newQs, newZs] = ordqz(RH, RR, Qs, Zs, clus_order);

newRH = newRH(1:new_nstp,1:new_nstp);
newRR = newRR(1:new_nstp,1:new_nstp);
new_Y = Y(:,1:nmax) * newQs(1:new_nstp,:)';
new_Z = Z * newZs(:,1:new_nstp);

tmp_vec          = newZs(end,1:new_nstp)';
[ tau, alpha, v] = householder_vector(tmp_vec(1,1),tmp_vec(2:new_nstp,1));
alpha            = conj(alpha);

% KZ    = mtx_K_prod_vec( mtx_G, mtx_M2, mtx_F, mtx_M1, new_Z(:,1:new_nstp), shift );
% Err_K = KZ - new_Y(:,1:new_nstp) * newRH(1:new_nstp,1:new_nstp) - ...
%         H(1+nmax,nmax) * Y(:,nmax+1) * tmp_vec';
% fprintf('Err_invK = %14.4e \n', norm(Err_K));
% NZ    = mtx_N_prod_vec( mtx_G, mtx_M2, mtx_F, mtx_M1, new_Z(:,1:new_nstp), shift );
% Err_N = NZ - new_Y(:,1:new_nstp) * newRR(1:new_nstp,1:new_nstp);
% fprintf('Err_invN = %14.4e \n', norm(Err_N));

if ( abs(alpha) > eps )
%
% Update new_Z by new_Z * ( I - tau * [1; v] * [ 1 v' ]
%
    vec_v = [ 1; v ];
    new_Z = new_Z - tau * (new_Z * vec_v) * vec_v';
    new_Z = [ new_Z(:,2:new_nstp) new_Z(:,1) ];

    newRH = newRH - tau * (newRH * vec_v) * vec_v';
    newRH = [ newRH(:,2:new_nstp) newRH(:,1) ];

    newRR = newRR - tau * (newRR * vec_v) * vec_v';
    newRR = [ newRR(:,2:new_nstp) newRR(:,1) ];

%     KZ    = mtx_K_prod_vec( mtx_G, mtx_M2, mtx_F, mtx_M1, new_Z(:,1:new_nstp), shift );
%     Err_K = KZ - new_Y(:,1:new_nstp) * newRH(1:new_nstp,1:new_nstp) - ...
%         H(1+nmax,nmax) * Y(:,nmax+1) * [zeros(1,new_nstp-1) 1];
%     fprintf('Err_invK = %14.4e \n', norm(Err_K));
%     NZ    = mtx_N_prod_vec( mtx_G, mtx_M2, mtx_F, mtx_M1, new_Z(:,1:new_nstp), shift );
%     Err_N = NZ - new_Y(:,1:new_nstp) * newRR(1:new_nstp,1:new_nstp);
%     fprintf('Err_invN = %14.4e \n', norm(Err_N));

    [new_H, new_R, Z_new, Y_new] = Schur_hess_tri_form( newRH, newRR, new_Z, new_Y );
else
    new_H = newRH;
    new_R = newRR;
    Z_new = new_Z;
    Y_new = new_Y;
end