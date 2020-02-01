function vec_y = FAME_Matrix_Vector_Production_Isotropic_shift_invAr_Simple( vec_x, B, Nx, Ny, Nz, N, Sigma_r, Pi_Qr, Pi_Qrs, D_k, D_ks, LS_tol, shift )
    global inner_iter inner_count inner_cpu_time 
    inner_count = inner_count + 1;
 
    vec_y = Sigma_r.\vec_x;
    time_start_pcg = tic;
    
%     [ vec_y, flag, res, inner_iter(inner_count), resvec] = bicgstab(  @(x) FAME_Matrix_Vector_Production_Isotropic_Ar_Simple(x, B, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_k, D_ks) - shift * ((Sigma_r.^2).\ x),...
%    vec_y, LS_tol, 100000 );
    [ vec_y, flag, res, inner_iter(inner_count), resvec] = minres(  @(x) FAME_Matrix_Vector_Production_Isotropic_Ar_Simple(x, B, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_k, D_ks) - shift * ((Sigma_r.^2).\ x),...
   vec_y, LS_tol, 10000 );
%     [ vec_y, flag, res, inner_iter(inner_count), resvec] = minres(  @(x) FAME_Matrix_Vector_Production_Isotropic_Ar_Simple_test(x, B, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_k, D_ks, shift) - (Sigma_r.^2).\ x ,...
%    vec_y, LS_tol, 10000);
    
%     [ vec_y, flag, res, inner_iter(inner_count), resvec] = minres(  @(x) FAME_Matrix_Vector_Production_Isotropic_Ar_test(x, B, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_k, D_ks, Sigma_r, shift),...
%    vec_y, LS_tol, 10000);
res
%     inner_iter(inner_count)
    inner_cpu_time(inner_count) = toc(time_start_pcg);
    
    vec_y = Sigma_r.\vec_y;
end


function vec_y = FAME_Matrix_Vector_Production_Isotropic_Ar_Simple_test(vec_x, B, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_k, D_ks, shift)
    vec_y = FAME_Matrix_Vector_Production_Qr_Simple(vec_x, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_k, D_ks, 'normal');
    vec_y = (1/shift) * vec_y;
%     vec_y = invB_eps.*vec_y;
    vec_y = FAME_Matrix_Vector_Production_invB_Isotropic(vec_y, B);

    vec_y = FAME_Matrix_Vector_Production_Qr_Simple(vec_y, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_k, D_ks, 'hermitian');
end