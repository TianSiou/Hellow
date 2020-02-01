  function vec_y = FAME_Matrix_Vector_Production_Isotropic_shift_invAr_General( vec_x, B, Nx, Ny, Nz, N, Sigma_r, Pi_Qr, Pi_Qrs, D_kx, D_ky, D_kz, F_ky, F_kz, LS_tol, shift)
    global inner_iter inner_count inner_cpu_time res
    inner_count = inner_count + 1;

    vec_y = Sigma_r.\vec_x;
    tic
%     time_start_pcg = tic;

%     [ vec_y, flag, res, inner_iter(inner_count), resvec] = bicgstab(  @(x) FAME_Matrix_Vector_Production_Isotropic_Ar_General(x, B, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_kx, D_ky, D_kz, F_ky, F_kz) - shift * ((Sigma_r.^2).\ x),...
%    vec_y, LS_tol, 100000 );
    [ vec_y, flag, res, inner_iter(inner_count), resvec] = minres(  @(x) FAME_Matrix_Vector_Production_Isotropic_Ar_General(x, B, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_kx, D_ky, D_kz, F_ky, F_kz) - shift * ((Sigma_r.^2).\ x),...
   vec_y, LS_tol, 10000 );

%     [ vec_y, flag, res, inner_iter(inner_count), resvec] = minres(  @(x) FAME_Matrix_Vector_Production_Isotropic_Ar_General(x, B, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_kx, D_ky, D_kz, F_ky, F_kz) - shift * ((Sigma_r.^2).\ x),...
%    vec_y, LS_tol, 10000, @(x) FAME_Matrix_Vector_Production_Isotropic_precondition_General( x, B, Nx, Ny, Nz, N, Sigma_r, Pi_Qr, Pi_Qrs, D_kx, D_ky, D_kz, F_ky, F_kz, 1e-10, shift));
  
res
    inner_iter(inner_count) 
%     inner_cpu_time(inner_count) = toc(time_start_pcg);
inner_cpu_time(inner_count) = toc
    vec_y = Sigma_r.\vec_y;
  end
  
  function vec_y = FAME_Matrix_Vector_Production_Isotropic_precondition_General( vec_x, B, Nx, Ny, Nz, N, Sigma_r, Pi_Qr, Pi_Qrs, D_kx, D_ky, D_kz, F_ky, F_kz, LS_tol, shift)
%     global  innerCG_count innerCG_cpu_time 
%     innerCG_count = innerCG_count + 1;

    [ vec_y, ~, ~, ~ ] = pcg(  @(x) FAME_Matrix_Vector_Production_Isotropic_Ar_General(x, B, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_kx, D_ky, D_kz, F_ky, F_kz) + shift * ((Sigma_r.^2).\ x),...
              vec_x, LS_tol, 1000 );
%     tic
%     [ vec_y, ~, ~, ~ ] = pcg(  @(x) FAME_Matrix_Vector_Production_Isotropic_Ar_General(x, B, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_kx, D_ky, D_kz, F_ky, F_kz),...
%               vec_x, LS_tol, 1000 );
%     innerCG_cpu_time(innerCG_count) = toc

end