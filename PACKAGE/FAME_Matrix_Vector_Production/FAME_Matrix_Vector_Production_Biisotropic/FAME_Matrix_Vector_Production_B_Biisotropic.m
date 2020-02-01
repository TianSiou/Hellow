function vec_y = FAME_Matrix_Vector_Production_B_Biisotropic(vec_x, B)
    vec_x_ele = vec_x(1:end/2);
    vec_x_mag = vec_x(end/2+1:end);
    
    vec_y = 1i*[B.B_zeta.*vec_x_ele + B.B_mu.*vec_x_mag; -B.B_eps.*vec_x_ele - B.B_xi.*vec_x_mag ];

end