function vec_y = FAME_Matrix_Vector_Production_Biisotropic_Ar_General_fem(vec_x, B, Nx, Ny, Nz, N, Pi_Qr, Pi_Pr, Pi_Qrs, Pi_Prs, D_kx, D_ky, D_kz)
    
    vec_x_ele = vec_x(1:end/2);
    vec_x_mag = vec_x(end/2+1:end);

    vec_y_ele = FAME_Matrix_Vector_Production_Pr_General(vec_x_ele, Nx, Ny, Nz, N, Pi_Pr, Pi_Prs, D_kx, D_ky, D_kz, 'normal');
    vec_y_mag = FAME_Matrix_Vector_Production_Qr_General(vec_x_mag, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_kx, D_ky, D_kz, 'normal');
    
    temp_vec_y_ele = B.Xi*vec_y_ele + vec_y_mag;
    temp_vec_y_mag = -vec_y_ele;
    
    vec_y_ele = FAME_Matrix_Vector_Production_invPhi_Biisotropic_fem(temp_vec_y_ele, B);
%     temp_vec_y_ele_1 = temp_vec_y_ele(1:end/3);
%     temp_vec_y_ele_2 = temp_vec_y_ele(end/3+1:2*end/3);
%     temp_vec_y_ele_3 = temp_vec_y_ele(2*end/3+1:end);
%     if strcmp(B.lssvr, 'mldivide')
%         vec_y_ele = B.Phi\temp_vec_y_ele;
%     elseif strcmp(B.lssvr, 'pcg')
%         time_start_pcg_B = tic;
%         [ vec_y_ele, ~, ~, inner_inner_iter(inner_inner_count) ] = pcg( B.Phi, temp_vec_y_ele, 1e-12, 1000);
%         inner_inner_cpu_time(inner_inner_count) = toc(time_start_pcg_B);
%     elseif isfield(B,'P_1')
%         temp_vec_y_ele_1 = B.L_1\temp_vec_y_ele_1(B.P_1);
%         temp_vec_y_ele_1 = B.U_1\temp_vec_y_ele_1;
%         vec_y_ele_1      = temp_vec_y_ele_1(B.invP_1);
%         temp_vec_y_ele_2 = B.L_2\temp_vec_y_ele_2(B.P_2);
%         temp_vec_y_ele_2 = B.U_2\temp_vec_y_ele_2;
%         vec_y_ele_2      = temp_vec_y_ele_2(B.invP_2);
%         temp_vec_y_ele_3 = B.L_3\temp_vec_y_ele_3(B.P_3);
%         temp_vec_y_ele_3 = B.U_3\temp_vec_y_ele_3;
%         vec_y_ele_3      = temp_vec_y_ele_3(B.invP_3);
%         
%         vec_y_ele = [vec_y_ele_1;vec_y_ele_2;vec_y_ele_3];
%     else
%         temp_vec_y_ele_1 = B.L_1\temp_vec_y_ele_1;
%         vec_y_ele_1      = B.U_1\temp_vec_y_ele_1;
%         temp_vec_y_ele_2 = B.L_2\temp_vec_y_ele_2;
%         vec_y_ele_2      = B.U_2\temp_vec_y_ele_2;
%         temp_vec_y_ele_3 = B.L_3\temp_vec_y_ele_3;
%         vec_y_ele_3      = B.U_3\temp_vec_y_ele_3;
%         
%         vec_y_ele = [vec_y_ele_1;vec_y_ele_2;vec_y_ele_3];
%     end

    vec_y_mag = temp_vec_y_mag;
    
    temp_vec_y_ele = B.Zeta*vec_y_ele - vec_y_mag;
    temp_vec_y_mag = vec_y_ele;
    
    vec_y_ele = FAME_Matrix_Vector_Production_Pr_General(temp_vec_y_ele, Nx, Ny, Nz, N, Pi_Pr, Pi_Prs, D_kx, D_ky, D_kz, 'hermitian');
    vec_y_mag = FAME_Matrix_Vector_Production_Qr_General(temp_vec_y_mag, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_kx, D_ky, D_kz, 'hermitian');    
    
    vec_y = [vec_y_ele; vec_y_mag];
end