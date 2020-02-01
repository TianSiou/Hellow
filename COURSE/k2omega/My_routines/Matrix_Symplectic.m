function mtx_symplecitc  = Matrix_Symplectic( mtx_palindromic, N )
    O_3N = sparse(3*N,3*N);
    I_3N = speye(3*N);
    
    mtx_symplecitc.fun_M = @(freq) [mtx_palindromic.fun_A(freq), O_3N; mtx_palindromic.fun_Q(freq), -I_3N];
    mtx_symplecitc.fun_L = @(freq) [O_3N, I_3N; mtx_palindromic.fun_A(freq).', O_3N];
end