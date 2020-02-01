function vec_y = FAME_Matrix_Vector_Production_invB_Isotropic(vec_x, B)
    vec_y = B.invB_eps.*vec_x;
end