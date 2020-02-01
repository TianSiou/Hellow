function flag = check_correction(N, mtx, fun_Tp, fun_Tps, option)

    vec_x = rand(6*N,1);
    % 測試 invN_hat 的 function handle 是否正確
    vec_y = Matrix_Vector_Production_invN_hat( vec_x, mtx.mtx_gyroscopic, mtx.mtx_palindromic, fun_Tp, fun_Tps, option);
    test_invN_hat = norm( vec_x - mtx.mtx_SH.N_hat(option.shift_palindromic)*vec_y )
    fprintf('\nCorrection of ''Matrix_Vector_Production_invN_hat'': ')
    if test_invN_hat < 1e-8
        fprintf('PASS\n');
        flag(1) = 0;
    else
        fprintf('FAIL\n');
        flag(1) = 1;
    end
    % 測試 invN2_hat 的 function handle 是否正確
    vec_y = Matrix_Vector_Production_invN2_hat(vec_x, mtx.mtx_palindromic, fun_Tp, fun_Tps, option);
    test_invN2_hat = norm( vec_x - mtx.mtx_SH.N2_hat(option.shift_palindromic)*vec_y )
    fprintf('\nCorrection of ''Matrix_Vector_Production_invN2_hat'': ')
    if test_invN2_hat < 1e-8
        fprintf('PASS\n');
        flag(1) = 0;
    else
        fprintf('FAIL\n');
        flag(1) = 1;
    end
    
    % 測試 W_hat 的 function handle 是否正確
    vec_y = mtx.mtx_SH.N2_hat(option.shift_palindromic)*vec_x;
    vec_y = Matrix_Vector_Production_W_hat(vec_y, mtx.mtx_palindromic, fun_Tp, fun_Tps, option);
    vec_y = mtx.mtx_SH.N1_hat(option.shift_palindromic)*vec_y;
    test_W_hat = norm( vec_y -  mtx.mtx_SH.K_hat(option.shift_palindromic)*vec_x)
    fprintf('\nCorrection of ''Matrix_Vector_Production_W_hat'': ')
    if test_W_hat < 1e-8
        fprintf('PASS\n');
        flag(1) = 0;
    else
        fprintf('FAIL\n');
        flag(1) = 1;
    end
end