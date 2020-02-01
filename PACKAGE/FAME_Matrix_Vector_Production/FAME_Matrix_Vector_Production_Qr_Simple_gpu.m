function dvec_y = FAME_Matrix_Vector_Production_Qr_Simple_gpu(dvec_x, Nx, Ny, Nz, N, dPi_Qr, dPi_Qrs, dD_k, dD_ks, mode)
    if strcmp(mode,'normal') == 1
        dvec_y = dPi_Qr*dvec_x;
        dvec_y = FAME_Matrix_Vector_Production_IFFT_Triple_Simple_gpu(dvec_y,Nx,Ny,Nz,N,dD_k);
    elseif strcmp(mode,'hermitian') == 1
        dvec_y = FAME_Matrix_Vector_Production_FFT_Triple_Simple_gpu(dvec_x,Nx,Ny,Nz,N,dD_ks);
        dvec_y = dPi_Qrs*dvec_y;
    end
end