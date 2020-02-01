function vec_y = Matrix_Vector_Production_invN2_hat(vec_x, mtx_palindromic, fun_invPp, option)
% reture vec_y = inv(N_hat)*vec_x
    nu_0 = option.shift_palindromic;

    vec_y1 =  vec_x(1:end/2) + (mtx_palindromic.Q - nu_0*mtx_palindromic.A)*vec_x(end/2+1:end);
    vec_y2 =  -vec_x(end/2+1:end);

    vec_y1 = Matrix_Vector_Production_invQp( vec_y1, 'transpose', mtx_palindromic, fun_invPp, option);
    
    vec_y2 = -nu_0*vec_y1 + vec_y2;
   
    vec_y = [ vec_y2;
             -vec_y1];
end