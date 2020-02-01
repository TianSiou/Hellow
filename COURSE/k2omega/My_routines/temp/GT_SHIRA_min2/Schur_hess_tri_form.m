function [mtx_H, mtx_R, Z_new, Y_new] = Schur_hess_tri_form( mtx_A, mtx_B, Z_old, Y_old )
% 
%        mtx_K * Z_old = Y_old * mtx_A + alpha * y * [ 0 ... 0 1 ]
%        mtx_N * Z_old = Y_old * mtx_B
%
% Produce an upper Hessenberg matrix mtx_H, an upper triangular matrix
% mtx_R, and unitary matrices mtx_Q and mtx_Z such that
%     mtx_Q' * mtx_A * mtx_Z = mtx_H,  mtx_Q' * mtx_B * mtx_Z = mtx_R
% and
%     [ 0 ... 0 1 ] * mtx_Z = [ 0 ... 0 1 ].
%
% Update
%     Z_new = Z_old * mtx_Z   and   Y_new = Y_old * mtx_Q
% so that
%        mtx_K * Z_new = Y_new * mtx_H + alpha * y * [ 0 ... 0 1 ]
%        mtx_N * Z_new = Y_new * mtx_R
%
n              = size(mtx_A,1);
mtx_Z          = eye(n);
[mtx_Q, mtx_R] = qr(mtx_B);
mtx_H          = mtx_Q' * mtx_A;
Y_new          = Y_old * mtx_Q;
Z_new          = Z_old;

for ii = n:-1:3
    for jj = 1:ii-2
        [G,Y]                 = planerot([mtx_H(ii,jj+1); mtx_H(ii,jj)]);
        mtx_H(ii,jj+1)        = Y(1);
        mtx_H(ii,jj)          = 0;
        mtx_G                 = [ G(2,2)  G(1,2); G(2,1)  G(1,1) ];
        mtx_H(1:ii-1,jj:jj+1) = mtx_H(1:ii-1,jj:jj+1) * mtx_G;
        mtx_R(1:jj+1,jj:jj+1) = mtx_R(1:jj+1,jj:jj+1) * mtx_G;
        mtx_Z(:,jj:jj+1)      = mtx_Z(:,jj:jj+1) * mtx_G;
        Z_new(:,jj:jj+1)      = Z_new(:,jj:jj+1) * mtx_G;
        %
        % Eliminate fill-in mtx_R(jj+1,jj) by Givens rotation from left 
        %
        [G,Y]                 = planerot(mtx_R(jj:jj+1,jj));
        mtx_R(jj,jj)          = Y(1);
        mtx_R(jj+1,jj)        = 0;
        mtx_R(jj:jj+1,jj+1:n) = G * mtx_R(jj:jj+1,jj+1:n);
        mtx_H(jj:jj+1,:)      = G * mtx_H(jj:jj+1,:);
        mtx_Q(:,jj:jj+1)      = mtx_Q(:,jj:jj+1) * G';
        Y_new(:,jj:jj+1)      = Y_new(:,jj:jj+1) * G';
    end
end
% Err_Z = Z_new - Z_old * mtx_Z;
% Err_Y = Y_new - Y_old * mtx_Q;
% fprintf('Error of Z = %12.4e, Error of Y = %12.4e\n',norm(Err_Z), norm(Err_Y));