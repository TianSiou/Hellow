function vec_y = FAME_Matrix_Vector_Production_B_Isotropic(vec_x, B)
    vec_y = B.B_eps.*vec_x;
end