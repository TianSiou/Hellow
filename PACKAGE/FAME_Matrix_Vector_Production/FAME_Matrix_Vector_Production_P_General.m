function vec_y = FAME_Matrix_Vector_Production_P_General(vec_x, Nx, Ny, Nz, N, Pi_P, Pi_Ps, D_kx, D_ky, D_kz, mode)
    if strcmp(mode,'normal') == 1
        vec_y = Pi_P*vec_x;
        vec_y = FAME_Matrix_Vector_Production_IFFT_Triple_General(vec_y,Nx,Ny,Nz,N, D_kx, D_ky, D_kz);
    elseif strcmp(mode,'hermitian') == 1
        vec_y = FAME_Matrix_Vector_Production_FFT_Triple_General(vec_x,Nx,Ny,Nz,N, D_kx, D_ky, D_kz);
        vec_y = Pi_Ps*vec_y;
    end
end