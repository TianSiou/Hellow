function vec_y = FAME_Matrix_Vector_Production_PML_Ahat_General(vec_x, freq, B_ele_tilde, B_mag_tilde, Avg_B_ele_tilde, Avg_B_mag_tilde, Normal_factor, Nx, Ny, Nz, N, Sigma_sqrt, Pi_Q, Pi_Qs, Pi_P, Pi_Ps, D_kx, D_ky, D_kz)
    vec_y = [ Avg_B_ele_tilde - B_ele_tilde ;
              Avg_B_mag_tilde - B_mag_tilde ].*vec_x;
    vec_y = 1i*freq*FAME_Matrix_Vector_Production_PML_invM_Simple(vec_y, freq, Avg_B_ele_tilde, Avg_B_mag_tilde, Normal_factor, Nx, Ny, Nz, N, Sigma_sqrt, Pi_Q, Pi_Qs, Pi_P, Pi_Ps, D_kx, D_ky, D_kz);
    vec_y = vec_x + vec_y;
end