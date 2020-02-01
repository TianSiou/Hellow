function funhand = Funhand_Generate(C, K_tilde, B, mtx, Par_mesh, Par_lattice, Lambdas, freq, option, N)
    funhand.gyroscopic.M  = @(vec_x) C'*(C*vec_x) - freq^2*(B.B_eps.*vec_x);
    funhand.gyroscopic.G  = @(vec_x) -K_tilde'*(C*vec_x) + C'*(K_tilde*vec_x);
    funhand.gyroscopic.K  = @(vec_x) -K_tilde'*(K_tilde*vec_x);
    funhand.palindromic.Qg = @(vec_x,mu_0) mu_0^2*funhand.gyroscopic.M(vec_x) + mu_0*funhand.gyroscopic.G(vec_x) + funhand.gyroscopic.K(vec_x);
    funhand.palindromic.A  = @(vec_x,transp_flag)  Matrix_Vector_Production_A(vec_x, transp_flag, funhand);
    funhand.palindromic.Q  = @(vec_x) 0.5*(funhand.gyroscopic.M(vec_x) - funhand.gyroscopic.K(vec_x));
    funhand.palindromic.Qp  = @(vec_x,nu_0) nu_0^2*funhand.palindromic.A(vec_x,'transp') - nu_0*funhand.palindromic.Q(vec_x) + funhand.palindromic.A(vec_x,'normal');
    funhand.palindromic.Qpt = @(vec_x,nu_0) nu_0^2*funhand.palindromic.A(vec_x,'normal') - nu_0*funhand.palindromic.Q(vec_x) + funhand.palindromic.A(vec_x,'transp');
    
%     Dp.d11 = mtx.mtx_palindromic.fun_Dp.fun_d11(option.shift_palindromic,freq,option.a_eps);
%     Dp.d12 = mtx.mtx_palindromic.fun_Dp.fun_d12(option.shift_palindromic);
%     Dp.d13 = mtx.mtx_palindromic.fun_Dp.fun_d13(option.shift_palindromic);
%     Dp.d21 = mtx.mtx_palindromic.fun_Dp.fun_d21(option.shift_palindromic);
%     Dp.d22 = mtx.mtx_palindromic.fun_Dp.fun_d22(option.shift_palindromic,freq,option.a_eps);
%     Dp.d23 = mtx.mtx_palindromic.fun_Dp.fun_d23(option.shift_palindromic);
%     Dp.d31 = mtx.mtx_palindromic.fun_Dp.fun_d31(option.shift_palindromic);
%     Dp.d32 = mtx.mtx_palindromic.fun_Dp.fun_d32(option.shift_palindromic);
%     Dp.d33 = mtx.mtx_palindromic.fun_Dp.fun_d33(option.shift_palindromic,freq,option.a_eps);
%     % 建構 invDp 的矩陣向量乘法的 function handle
%     funhand.invD_fun = Matrix_Inverse_BlockMatrix(Dp);
    % 建構 FFT-based 矩陣向量乘法 Tp 與 Tps 的 function handle
    switch Par_lattice.lattice_type
        case {'simple_cubic', 'primitive_orthorhombic', 'primitive_tetragonal'}
            funhand.fun_Tp  = @(vec_x) FAME_Matrix_Vector_Production_IFFT_Triple_Simple(vec_x, Par_mesh.grid_num(1), Par_mesh.grid_num(2), Par_mesh.grid_num(3), N, Lambdas.D_k);
            funhand.fun_Tps = @(vec_x) FAME_Matrix_Vector_Production_FFT_Triple_Simple(vec_x,Par_mesh.grid_num(1),Par_mesh.grid_num(2),Par_mesh.grid_num(3),N,Lambdas.D_ks);
        otherwise
            funhand.fun_Tp  = @(vec_x) FAME_Matrix_Vector_Production_IFFT_Triple_General(vec_x,Par_mesh.grid_num(1),Par_mesh.grid_num(2),Par_mesh.grid_num(3),N,Lambdas.D_kx,Lambdas.F_ky,Lambdas.F_kz);
            funhand.fun_Tps = @(vec_x) FAME_Matrix_Vector_Production_FFT_Triple_General(vec_x,Par_mesh.grid_num(1),Par_mesh.grid_num(2),Par_mesh.grid_num(3),N,Lambdas.D_kx,Lambdas.D_ky,Lambdas.D_kz);
    end  
    
    % 建構 W_hat 的 function handle
%     funhand.fun_W_hat = @(vec_x) Matrix_Vector_Production_W_hat(vec_x, mtx.mtx_palindromic, funhand.fun_Tp, funhand.fun_Tps, option);
    % 建構 invN2_hat 的 function handle
%     funhand.fun_invN2_hat = @(vec_x) Matrix_Vector_Production_invN2_hat(vec_x, mtx.mtx_palindromic, funhand.fun_Tp, funhand.fun_Tps, option);
    % 建構 invPp 的 function handle
    funhand.fun_invPp = @(vec_x,mode) Matrix_Vector_Production_invPp(vec_x, mode, mtx.mtx_palindromic, option.linear_system.alpha, funhand.fun_Tp, funhand.fun_Tps);
%     funhand.fun_Pp    = @(vec_x) Matrix_Vector_Production_Pp(vec_x, mtx.mtx_palindromic.Sigma_pp, mtx.mtx_palindromic.Lambda_c, funhand.fun_Tp, funhand.fun_Tps);
    % 建構 invQp 的 function handle
%     funhand.fun_invQp = @(vec_x,mode) Matrix_Vector_Production_invQp( vec_x, mode, mtx.mtx_palindromic, funhand.fun_invPp, option);
    % 建構 invN_hat_K_hat 的 function handle
%     funhand.fun_invN_hat_K_hat = @(vec_x) Matrix_Vector_Production_invN_hat( mtx.mtx_SH.K_hat*vec_x, ...
%                                                                              mtx.mtx_gyroscopic, mtx.mtx_palindromic, funhand.fun_invPp, option);      
    
                                                                         
%     vec_x = rand(3*N,1);
%     vec_y = funhand.fun_invPp(vec_x, 'normal');
%     vec_y = funhand.fun_invPp(vec_x, 'normal');
%     test_Pp  = norm(vec_x - funhand.fun_Pp(vec_y))
%     test_Pp  = norm(vec_x - funhand.fun_Pp(vec_y))
%     test_Ppt = norm(vec_x - mtx.mtx_palindromic.Ppt_funhand(funhand.fun_invPp(vec_x, 'trans'), option.shift_palindromic,option.alpha))
%     test = norm(vec_x - mtx.mtx_palindromic.PpT*funhand.fun_invPp(vec_x, 'trans'))
%     test = norm(vec_x - mtx.mtx_palindromic.Qp*funhand.fun_invQp(vec_x, 'normal'))
%     test = norm(vec_x - mtx.mtx_palindromic.QpT*funhand.fun_invQp(vec_x, 'trans'))
    % 檢驗 function handle 正確性
%     check_correction(N,mtx, funhand.fun_Tp, funhand.fun_Tps, option);
%     funhand.fun_test = @(vec_x,mode) Matrix_Vector_Production_test( vec_x, mode, mtx.mtx_palindromic, funhand.fun_invPp, option);
%     test = norm(vec_x - funhand.fun_test(vec_x, 'normal'))
end

function vec_y = Matrix_Vector_Production_Pp(vec_x, mtx_Sigma_pp, mtx_Lambda_c, fun_Tp, fun_Tps)
    vec_y = fun_Tps(vec_x);
    vec_y = vec_y.*mtx_Sigma_pp - mtx_Lambda_c*(mtx_Lambda_c'*vec_y);
    vec_y = fun_Tp(vec_y);
end

function vec_y = Matrix_Vector_Production_A(vec_x, transp_flag, funhand)
    if strcmp(transp_flag,'transp')
        vec_y  = 0.25*(funhand.gyroscopic.M(vec_x) + funhand.gyroscopic.G(vec_x) + funhand.gyroscopic.K(vec_x));
    else
        vec_y  = 0.25*(funhand.gyroscopic.M(vec_x) - funhand.gyroscopic.G(vec_x) + funhand.gyroscopic.K(vec_x));
    end
end