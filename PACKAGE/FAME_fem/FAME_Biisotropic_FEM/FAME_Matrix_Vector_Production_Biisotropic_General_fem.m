function vec_y = FAME_Matrix_Vector_Production_Biisotropic_General_fem( vec_x, B, Nx, Ny, Nz, N, inv_Sigma_r, Pi_Qr, Pi_Pr, Pi_Qrs, Pi_Prs, D_kx, D_ky, D_kz, LS_tol)
    vec_x_ele = vec_x(1:end/2);
    vec_x_mag = vec_x(end/2+1:end);

    vec_y_ele =  inv_Sigma_r.*vec_x_mag;
    vec_y_mag = -inv_Sigma_r.*vec_x_ele;
    
    vec_y = 1i*[vec_y_ele;vec_y_mag];
    
    vec_y = FAME_Matrix_Vector_Production_Biisotropic_invAr_General_fem( vec_y, B, Nx, Ny, Nz, N, Pi_Qr, Pi_Pr, Pi_Qrs, Pi_Prs, D_kx, D_ky, D_kz, LS_tol);
end