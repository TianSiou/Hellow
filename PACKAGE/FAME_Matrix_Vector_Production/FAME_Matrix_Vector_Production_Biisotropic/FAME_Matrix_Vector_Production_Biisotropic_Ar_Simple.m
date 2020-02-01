function vec_y = FAME_Matrix_Vector_Production_Biisotropic_Ar_Simple(vec_x, B, Nx, Ny, Nz, N, Pi_Qr, Pi_Pr, Pi_Qrs, Pi_Prs, D_k, D_ks)
    vec_x_ele = vec_x(1:end/2);
    vec_x_mag = vec_x(end/2+1:end);

    vec_y_ele = FAME_Matrix_Vector_Production_Pr_Simple(vec_x_ele, Nx, Ny, Nz, N, Pi_Pr, Pi_Prs, D_k, D_ks, 'normal');
    vec_y_mag = FAME_Matrix_Vector_Production_Qr_Simple(vec_x_mag, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_k, D_ks, 'normal');
    
    temp_vec_y_ele = B.B_zeta_s.*vec_y_ele + vec_y_mag;
    temp_vec_y_mag = -vec_y_ele;
    
%     vec_y_ele = invPhi.*temp_vec_y_ele;
    vec_y_ele = FAME_Matrix_Vector_Production_invPhi_Biisotropic(temp_vec_y_ele, B);
    vec_y_mag = temp_vec_y_mag;
    
    temp_vec_y_ele = B.B_zeta.*vec_y_ele - vec_y_mag;
    temp_vec_y_mag = vec_y_ele;
    
    vec_y_ele = FAME_Matrix_Vector_Production_Pr_Simple(temp_vec_y_ele, Nx, Ny, Nz, N, Pi_Pr, Pi_Prs, D_k, D_ks, 'hermitian');
    vec_y_mag = FAME_Matrix_Vector_Production_Qr_Simple(temp_vec_y_mag, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_k, D_ks, 'hermitian');    
    
    vec_y = [vec_y_ele; vec_y_mag];
end