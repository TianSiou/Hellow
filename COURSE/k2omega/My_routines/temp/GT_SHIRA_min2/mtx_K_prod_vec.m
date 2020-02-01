function MtxVec = mtx_K_prod_vec(mtx_A, vec, shift)
%
%                    | A    0   |
%   Compute lambda * |          | * x.
%                    | 0    A^T |
%  
if isa(mtx_A,'function_handle')
    n      = size(vec, 1) / 2;
    tmp1   = mtx_A(vec(1:n, :), 'notransp');
    tmp2   = mtx_A(vec((n+1):(2*n), :), 'transp');
    MtxVec = shift * [tmp1; tmp2];
    
else
    % Check matrix and right hand side vector inputs have appropriate sizes
    n      = size(mtx_A, 1);
    tmp1   = mtx_A * vec(1:n, :);
    tmp2   = (mtx_A.') * vec((n+1):(2*n), :);
    MtxVec = shift * [tmp1; tmp2];
    
end

end
