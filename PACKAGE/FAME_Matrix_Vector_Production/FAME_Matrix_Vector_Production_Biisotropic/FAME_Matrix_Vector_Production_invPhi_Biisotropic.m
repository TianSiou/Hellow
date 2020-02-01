function vec_y = FAME_Matrix_Vector_Production_invPhi_Biisotropic(vec_x, B)
    vec_y = B.invPhi.*vec_x;
end