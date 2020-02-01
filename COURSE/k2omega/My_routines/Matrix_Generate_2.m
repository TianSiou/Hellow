function mtx = Matrix_Generate_2(mtx, freq, nu_0, alpha)
% Generate frequency-dependent matrices
    mtx.mtx_gyroscopic.M    = mtx.mtx_gyroscopic.fun_M(freq);
%     mtx.mtx_gyroscopic.Qg   = mtx.mtx_gyroscopic.fun_Qg(nu_0,freq);
%     mtx.mtx_gyroscopic.Pg   = mtx.mtx_gyroscopic.fun_Pg(nu_0,alpha);
    
    mtx.mtx_SHH.H           = mtx.mtx_SHH.fun_H(freq);
    mtx.mtx_SHH.S           = mtx.mtx_SHH.fun_S(freq);
%     mtx.mtx_SHH.S_1         = mtx.mtx_SHH.fun_S_1(freq);
%     mtx.mtx_SHH.S_2         = mtx.mtx_SHH.fun_S_2(freq);
%     
    mtx.mtx_palindromic.A             = mtx.mtx_palindromic.fun_A(freq);
    mtx.mtx_palindromic.Q             = mtx.mtx_palindromic.fun_Q(freq);
%     mtx.mtx_palindromic.Qp            = mtx.mtx_palindromic.fun_Qp(nu_0,freq);
%     if nu_0 == -1
%         mtx.mtx_palindromic.Qp = mtx.C.'*mtx.C - freq^2*(mtx.B_eps);
%     end
%     mtx.mtx_palindromic.QpT           = mtx.mtx_palindromic.Qp.';
%     mtx.mtx_palindromic.Pp            = mtx.mtx_palindromic.fun_Pp(nu_0,alpha);
%     mtx.mtx_palindromic.PpT           = mtx.mtx_palindromic.Pp.';
    mtx.mtx_palindromic.Lambda_N1     = mtx.mtx_palindromic.fun_Lambda_N1(nu_0);
    mtx.mtx_palindromic.Lambda_N2     = mtx.mtx_palindromic.fun_Lambda_N2(nu_0);
    mtx.mtx_palindromic.Lambda_N1bar  = mtx.mtx_palindromic.fun_Lambda_N1bar(nu_0);
    mtx.mtx_palindromic.Lambda_N2bar  = mtx.mtx_palindromic.fun_Lambda_N2bar(nu_0);
    mtx.mtx_palindromic.invSigma_p    = mtx.mtx_palindromic.fun_invSigma_p(nu_0,alpha);
    mtx.mtx_palindromic.invSigma_pt   = mtx.mtx_palindromic.fun_invSigma_pt(nu_0,alpha);
    
    mtx.mtx_symplecitc.M    = mtx.mtx_symplecitc.fun_M(freq);
    mtx.mtx_symplecitc.L    = mtx.mtx_symplecitc.fun_L(freq);
% 
    mtx.mtx_ss.Ms       = mtx.mtx_ss.fun_Ms(freq);
    mtx.mtx_ss.Ls       = mtx.mtx_ss.fun_Ls(freq);
%   
    mtx.mtx_SH.K            = mtx.mtx_SH.fun_K(freq);
    mtx.mtx_SH.N            = mtx.mtx_SH.fun_N(freq);
%     mtx.mtx_SH.K_hat        = mtx.mtx_SH.fun_K_hat(nu_0,freq);
%     mtx.mtx_SH.N_hat        = mtx.mtx_SH.fun_N_hat(nu_0,freq);
%     mtx.mtx_SH.N1_hat       = mtx.mtx_SH.fun_N1_hat(nu_0,freq);
%     mtx.mtx_SH.N2_hat       = mtx.mtx_SH.fun_N2_hat(nu_0,freq);