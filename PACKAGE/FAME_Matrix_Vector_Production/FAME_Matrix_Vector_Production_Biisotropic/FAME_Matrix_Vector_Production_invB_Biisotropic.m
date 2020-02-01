function vec_y = FAME_Matrix_Vector_Production_invB_Biisotropic(vec_x, B)
    vec_x_ele = vec_x(1:end/2);
    vec_x_mag = vec_x(end/2+1:end);
    
    vec_y_ele = (-B.B_xi .*vec_x_ele - B.B_mu  .*vec_x_mag).*B.invPhi;
    vec_y_mag = ( B.B_eps.*vec_x_ele + B.B_zeta.*vec_x_mag).*B.invPhi;
    
    vec_y = -1i*[vec_y_ele;vec_y_mag];
end