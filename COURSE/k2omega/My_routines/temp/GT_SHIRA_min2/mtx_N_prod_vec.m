function MtxVec = mtx_N_prod_vec(mtx_A, mtx_Q, vec, shift)
%
% Compute (M^T - lambda L^T) J^T x
%
n = size(vec, 1) / 2;

if isa(mtx_Q,'function_handle')
    tmp1 = mtx_Q(vec(1:n, :));
else
    tmp1 = mtx_Q * vec(1:n, :);
end

if isa(mtx_A,'function_handle')
    tmp2   = mtx_A(vec(1:n, :), 'notransp');
    tmp3   = mtx_A(vec((n+1):(2*n), :), 'transp');
else
    tmp2 = mtx_A * vec(1:n, :);
    tmp3 = mtx_A.' * vec(n+1:2*n, :);
end
%     tmp1 = Q * vec(1:n, :);
%     tmp2 = A * vec(1:n, :);
%     tmp3 = A.' * vec(n+1:2*n, :);
    tmp  = [tmp1 - shift*tmp2 - tmp3; -vec(1:n, :) + shift*vec(n+1:2*n, :)];
    
%
% Compute (M - lambda L) J x
%
if isa(mtx_Q,'function_handle')
    tmp2 = mtx_Q(tmp(n+1:2*n, :));
else
    tmp2 = mtx_Q * tmp(n+1:2*n, :);
end

if isa(mtx_A,'function_handle')
    tmp1   = mtx_A(tmp(n+1:2*n, :), 'notransp');
    tmp3   = mtx_A(tmp(n+1:2*n, :), 'transp');
else
    tmp1 = mtx_A * tmp(n+1:2*n, :);
    tmp3 = mtx_A.' * tmp(n+1:2*n, :);
end

%     tmp1 = A * tmp(n+1:2*n, :);
%     tmp2 = Q * tmp(n+1:2*n, :);
%     tmp3 = A.' * tmp(n+1:2*n, :);
    MtxVec = [shift*tmp(1:n, :) + tmp1; tmp(1:n, :) + tmp2 - shift*tmp3 ];
end