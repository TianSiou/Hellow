function vec_y = FAME_Matrix_Vector_Production_Pr_General(vec_x, Nx, Ny, Nz, N, Pi_Pr, Pi_Prs, D_kx, D_ky, D_kz, F_ky, F_kz, mode)
    if strcmp(mode,'normal') == 1
        vec_y = Pi_Pr*vec_x;
        vec_y = FAME_Matrix_Vector_Production_IFFT_Triple_General(vec_y,Nx,Ny,Nz,N, D_kx, F_ky, F_kz);
    elseif strcmp(mode,'hermitian') == 1
        vec_y = FAME_Matrix_Vector_Production_FFT_Triple_General(vec_x,Nx,Ny,Nz,N, D_kx, D_ky, D_kz);
        vec_y = Pi_Prs*vec_y;
    end
end