function vec_y = Matrix_Vector_Production_invN_hat(vec_x, mtx_gyroscopic, mtx_palindromic, fun_invPp, option)
% reture vec_y = inv(N_hat)*vec_x
global iter_fun_invN_hat
iter_fun_invN_hat = iter_fun_invN_hat + 1;

    nu_0 = option.shift_palindromic;

    vec_y1 =  vec_x(1:end/2) - nu_0*vec_x(end/2+1:end);
    vec_y2 =  vec_x(end/2+1:end);

    vec_y1 = Matrix_Vector_Production_invQp( vec_y1, 'normal', mtx_palindromic, fun_invPp, option);

    vec_y1_temp = 0.5*nu_0*(mtx_gyroscopic.G*vec_y1) + vec_y2;
    vec_y2      = vec_y1;

    vec_y1 = Matrix_Vector_Production_invQp( vec_y1_temp, 'transpose', mtx_palindromic, fun_invPp, option);

    vec_y = -[nu_0*vec_y1 + vec_y2;
              vec_y1              ];
end