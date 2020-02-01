function vec_y = Matrix_Vector_Production_W_hat(vec_x, mtx_palindromic, fun_invPp, option)
% reture vec_y = W_hat*vec_x
global iter_fun_W_hat
iter_fun_W_hat = iter_fun_W_hat + 1;

    nu_0 = option.shift_palindromic;

    vec_y1 =   vec_x(1:end/2) + (mtx_palindromic.Q - nu_0*mtx_palindromic.A)*vec_x(end/2+1:end);
    vec_y2 =  -vec_x(end/2+1:end);

    vec_y1 = Matrix_Vector_Production_invQp( vec_y1, 'transpose', mtx_palindromic, fun_invPp, option);
    
    vec_y2 = -nu_0*vec_y1 + vec_y2;
    
    vec_y1_temp =  nu_0*mtx_palindromic.A  *vec_y2;
    vec_y2      = -nu_0*mtx_palindromic.A.'*vec_y1;
    
    vec_y1 = vec_y1_temp - nu_0*vec_y2;
   
    vec_y1 = Matrix_Vector_Production_invQp( vec_y1, 'normal', mtx_palindromic, fun_invPp, option);

    vec_y1_temp = (-mtx_palindromic.Q + nu_0*mtx_palindromic.A.')*vec_y1 + vec_y2;
    vec_y2 = vec_y1;
    
    vec_y = [vec_y1_temp;
             vec_y2    ];
end