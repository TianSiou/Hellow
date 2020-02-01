function vec_y = FAME_Matrix_Vector_Production_Anisotropic_Ar_Simple(vec_x, B, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_k, D_ks)
    vec_y = FAME_Matrix_Vector_Production_Qr_Simple(vec_x, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_k, D_ks, 'normal');
  
%     vec_y = invB_eps*vec_y;
%     vec_y = [ invB_eps.d11.*vec_y(1:N) + invB_eps.d12.*vec_y(N+1:2*N) + invB_eps.d13.*vec_y(2*N+1:3*N);
%               invB_eps.d21.*vec_y(1:N) + invB_eps.d22.*vec_y(N+1:2*N) + invB_eps.d23.*vec_y(2*N+1:3*N);
%               invB_eps.d31.*vec_y(1:N) + invB_eps.d32.*vec_y(N+1:2*N) + invB_eps.d33.*vec_y(2*N+1:3*N) ];
    vec_y = FAME_Matrix_Vector_Production_invB_Anisotropic(vec_y, B, N);

    vec_y = FAME_Matrix_Vector_Production_Qr_Simple(vec_y, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_k, D_ks, 'hermitian');
end