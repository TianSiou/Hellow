function vec_y = FAME_Matrix_Vector_Production_Isotropic_Ar_Simple(vec_x, B, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_k, D_ks)
    vec_y = FAME_Matrix_Vector_Production_Qr_Simple(vec_x, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_k, D_ks, 'normal');
%     vec_y = (1/shift) * vec_y;
%     vec_y = invB_eps.*vec_y;
    vec_y = FAME_Matrix_Vector_Production_invB_Isotropic(vec_y, B);

    vec_y = FAME_Matrix_Vector_Production_Qr_Simple(vec_y, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_k, D_ks, 'hermitian');
end