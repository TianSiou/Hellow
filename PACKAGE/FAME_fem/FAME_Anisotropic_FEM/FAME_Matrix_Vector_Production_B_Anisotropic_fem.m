function vec_y = FAME_Matrix_Vector_Production_B_Anisotropic_fem(vec_x, B)
    vec_y = B.B*vec_x;
end