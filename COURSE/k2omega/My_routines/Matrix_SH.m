function mtx_SH = Matrix_SH( mtx_palindromic, mtx_symplecitc, N )
    O_3N = sparse(3*N,3*N);
    I_3N = speye(3*N);
    J = [ O_3N,I_3N ;
         -I_3N,O_3N];
    
    mtx_SH.fun_K = @(freq) [-mtx_palindromic.fun_Q(freq)                ,  mtx_palindromic.fun_A(freq).'-mtx_palindromic.fun_A(freq) ; 
                             mtx_palindromic.fun_A(freq)-mtx_palindromic.fun_A(freq).', -mtx_palindromic.fun_Q(freq)                              ];
    mtx_SH.fun_N = @(freq) [-mtx_palindromic.fun_A(freq),                         O_3N   ; 
                                                O_3N, -mtx_palindromic.fun_A(freq).'];
    
    mtx_SH.fun_K_hat = @(sigma,freq) -sigma*mtx_SH.fun_N(freq);
    mtx_SH.fun_N_hat = @(sigma,freq) -sigma*(mtx_SH.fun_K(freq) - (sigma+1/sigma)*mtx_SH.fun_N(freq));
    mtx_SH.fun_N1_hat = @(sigma,freq) (mtx_symplecitc.fun_M(freq) - sigma*mtx_symplecitc.fun_L(freq))*J;
    mtx_SH.fun_N2_hat = @(sigma,freq) (mtx_symplecitc.fun_M(freq) - sigma*mtx_symplecitc.fun_L(freq)).'*J.';
    
%     mu_0 = rand(1) + rand(1)*1i;
%     test_N1N2 = norm( mtx_SH.N_hat(mu_0) - mtx_SH.N1_hat(mu_0)*mtx_SH.N2_hat(mu_0), 'fro' )
end