function vec_y = FAME_Matrix_Vector_Production_PML_M_General(vec_x, freq, Avg_B_ele_tilde, Avg_B_mag_tilde, Nx, Ny, Nz, N, Sigma_sqrt, Pi_Q, Pi_Qs, Pi_P, Pi_Ps, D_kx, D_ky, D_kz)
    % Preconditioner for matrix-vector production Ax

    vec_x1 = vec_x(1:end/2);
    vec_x2 = vec_x(end/2+1:end);
    
    vec_y1 = FAME_Matrix_Vector_Production_Q_General(vec_x1, Nx, Ny, Nz, N, Pi_Q, Pi_Qs, D_kx, D_ky, D_kz, 'hermitian');
    vec_y2 = FAME_Matrix_Vector_Production_P_General(vec_x2, Nx, Ny, Nz, N, Pi_P, Pi_Ps, D_kx, D_ky, D_kz, 'hermitian');
    
    y1_new = -1i*freq*Avg_B_ele_tilde*vec_y1 - Sigma_sqrt.*vec_y2;
    y2_new =  Sigma_sqrt.*vec_y1 - 1i*freq*Avg_B_mag_tilde*vec_y2;

    vec_y1 = FAME_Matrix_Vector_Production_Q_General(y1_new, Nx, Ny, Nz, N, Pi_Q, Pi_Qs, D_kx, D_ky, D_kz, 'normal');
    vec_y2 = FAME_Matrix_Vector_Production_P_General(y2_new, Nx, Ny, Nz, N, Pi_P, Pi_Ps, D_kx, D_ky, D_kz, 'normal');
        
    vec_y = [ vec_y1; vec_y2 ];
end