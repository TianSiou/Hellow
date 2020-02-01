function mtx_gyroscopic = Matrix_Gyroscopic( C, K_tilde, B )
    N = length(K_tilde);

    mtx_gyroscopic.K     = -K_tilde.'*K_tilde;
    mtx_gyroscopic.G     = C.'*K_tilde - K_tilde.'*C;
  
    mtx_gyroscopic.fun_M = @(freq) C.'*C - (freq^2)*spdiags(B.B_eps,0,N,N);
%     a_eps = sum(B.B_eps)/length(B.B_eps);
    
    mtx_gyroscopic.fun_Qg = @(sigma,freq)  sigma^2*mtx_gyroscopic.fun_M(freq)        + sigma*mtx_gyroscopic.G + mtx_gyroscopic.K;
    mtx_gyroscopic.fun_Pg = @(sigma,alpha) sigma^2*(C.'*C - alpha*speye(length(mtx_gyroscopic.K))) + sigma*mtx_gyroscopic.G + mtx_gyroscopic.K;
    
    if 0
        test = rand(1);
        norm(mtx_gyroscopic.fun_M(test) - mtx_gyroscopic.fun_M(test).','fro')
        norm(mtx_gyroscopic.G + mtx_gyroscopic.G.','fro')
        norm(mtx_gyroscopic.K - mtx_gyroscopic.K.','fro')
    end
end