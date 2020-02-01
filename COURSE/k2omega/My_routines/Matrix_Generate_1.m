function mtx = Matrix_Generate_1(C, C_1, C_2, C_3, K_tilde, K_1_tilde, K_2_tilde, K_3_tilde, B, Lambdas, N)

% mtx.mtx_palindromic = Matrix_Palindromic( C, C_1, C_2, C_3, K_tilde, K_1_tilde, K_2_tilde, K_3_tilde, [], Lambdas );    

%     mtx.C   = C;    mtx.K_tilde   = K_tilde;
%     mtx.C_1 = C_1;  mtx.K_1_tilde = K_tilde;
%     mtx.C_2 = C_2;  mtx.K_2_tilde = K_tilde;
%     mtx.C_3 = C_3;  mtx.K_3_tilde = K_tilde;
%     mtx.B_eps = spdiags(B.B_eps,0,3*N,3*N);

    mtx.mtx_gyroscopic  = Matrix_Gyroscopic( C, K_tilde, B );    
    mtx.mtx_SHH         = Matrix_SHH( mtx.mtx_gyroscopic);
    mtx.mtx_palindromic = Matrix_Palindromic( C, C_1, C_2, C_3, K_tilde, K_1_tilde, K_2_tilde, K_3_tilde, mtx.mtx_gyroscopic, Lambdas );    
    mtx.mtx_symplecitc  = Matrix_Symplectic( mtx.mtx_palindromic, N );    
    mtx.mtx_ss          = Matrix_SS( mtx.mtx_symplecitc, N ); 
    mtx.mtx_SH          = Matrix_SH( mtx.mtx_palindromic, mtx.mtx_symplecitc, N ); 
    
    % 建構 3x3 diagonal block matrix: invDp 和他的共扼轉置
%     mtx.mtx_palindromic.invDp  = mtx.mtx_palindromic.fun_invDp(option.shift_palindromic,option.a_eps );
%     mtx.mtx_palindromic.invDpT = mtx.mtx_palindromic.fun_invDpT(option.shift_palindromic,option.a_eps);
%     mtx.mtx_palindromic.Qp     = mtx.mtx_palindromic.fun_Qp(option.shift_palindromic);
%     mtx.mtx_palindromic.QpT    = mtx.mtx_palindromic.Qp.';
    
%     D_fun(1:N,    1:  N)          = spdiags(D.d11,0,N,N);
%     D_fun(1:N,  N+1:2*N)          = spdiags(D.d12,0,N,N);
%     D_fun(1:N,2*N+1:3*N)          = spdiags(D.d13,0,N,N);
%     D_fun(N+1:2*N,    1:  N)      = spdiags(D.d21,0,N,N);
%     D_fun(N+1:2*N,  N+1:2*N)      = spdiags(D.d22,0,N,N);
%     D_fun(N+1:2*N,2*N+1:3*N)      = spdiags(D.d23,0,N,N);
%     D_fun(2*N+1:3*N,    1:  N)    = spdiags(D.d31,0,N,N);
%     D_fun(2*N+1:3*N,  N+1:2*N)    = spdiags(D.d32,0,N,N);
%     D_fun(2*N+1:3*N,2*N+1:3*N)    = spdiags(D.d33,0,N,N);
    
%     vec_x = rand(3*N,1);
%     tic
%     vec_y = invD_fun(vec_x);
%     time_invD_fun = toc
%     tic
%     vec_y = D_fun\vec_x;
%     time_ldivide = toc
%     test_invDp = norm(vec_x - D_fun*vec_y )
%     test_invDp = norm(speye(3*N) - mtx.mtx_palindromic.invDp*mtx.mtx_palindromic.fun_Dp(option.shift_palindromic,option.a_eps),'fro')
end