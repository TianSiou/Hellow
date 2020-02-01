function vec_y = FAME_Matrix_Vector_Production_P_Simple(vec_x, Nx, Ny, Nz, N, Pi_P, Pi_Ps, D_k, D_ks, mode)
    if strcmp(mode,'normal') == 1
        vec_y = Pi_P*vec_x;
        vec_y = FAME_Matrix_Vector_Production_IFFT_Triple_Simple(vec_y,Nx,Ny,Nz,N,D_k);
    elseif strcmp(mode,'hermitian') == 1
        vec_y = FAME_Matrix_Vector_Production_FFT_Triple_Simple(vec_x,Nx,Ny,Nz,N,D_ks);
        vec_y = Pi_Ps*vec_y;
    end
end