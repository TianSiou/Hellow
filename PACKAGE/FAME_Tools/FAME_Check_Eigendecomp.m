function  [ flag_all, flag_C1, flag_C2, flag_C3, flag_CsC ] = FAME_Check_Eigendecomp(grid_num, lattice_type, C, Cs, C_1, C_2, C_3, Lambdas)
    flag_C1  = 1;
    flag_C2  = 1;
    flag_C3  = 1;
    flag_CsC = 1;
    flag_all = 1;
    switch lattice_type
        case {'simple_cubic', 'primitive_orthorhombic', 'primitive_tetragonal'}
            x    = rand(grid_num(1)*grid_num(2)*grid_num(3),1);
            Ts_x = FAME_Matrix_Vector_Production_FFT_Single_Simple(x,grid_num(1), grid_num(2), grid_num(3), grid_num(1)*grid_num(2)*grid_num(3),Lambdas.D_ks);
            
            lambdax_Ts_x = Lambdas.Lambda_x.*Ts_x;
            lambday_Ts_x = Lambdas.Lambda_y.*Ts_x;
            lambdaz_Ts_x = Lambdas.Lambda_z.*Ts_x;
            
            y_1 = FAME_Matrix_Vector_Production_IFFT_Single_Simple(lambdax_Ts_x,grid_num(1), grid_num(2), grid_num(3), grid_num(1)*grid_num(2)*grid_num(3),Lambdas.D_k);
            y_2 = FAME_Matrix_Vector_Production_IFFT_Single_Simple(lambday_Ts_x,grid_num(1), grid_num(2), grid_num(3), grid_num(1)*grid_num(2)*grid_num(3),Lambdas.D_k);
            y_3 = FAME_Matrix_Vector_Production_IFFT_Single_Simple(lambdaz_Ts_x,grid_num(1), grid_num(2), grid_num(3), grid_num(1)*grid_num(2)*grid_num(3),Lambdas.D_k);
            
            test_x = norm(C_1*x - y_1);
            test_y = norm(C_2*x - y_2);
            test_z = norm(C_3*x - y_3);
            
            x             = rand(3*grid_num(1)*grid_num(2)*grid_num(3),1);
            Qrs_x         = FAME_Matrix_Vector_Production_Qr_Simple(x, grid_num(1), grid_num(2), grid_num(3), grid_num(1)*grid_num(2)*grid_num(3), Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_k, Lambdas.D_ks, 'hermitian');
            Sigma_r_Qrs_x = (Lambdas.Sigma_r.^2).*Qrs_x;
            y             = FAME_Matrix_Vector_Production_Qr_Simple(Sigma_r_Qrs_x, grid_num(1), grid_num(2), grid_num(3), grid_num(1)*grid_num(2)*grid_num(3), Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_k, Lambdas.D_ks, 'normal');
                  
        otherwise  
            x    = rand(grid_num(1)*grid_num(2)*grid_num(3),1);
            Ts_x = FAME_Matrix_Vector_Production_FFT_Single_General(x,grid_num(1), grid_num(2), grid_num(3), grid_num(1)*grid_num(2)*grid_num(3),Lambdas.D_kx,Lambdas.D_ky,Lambdas.D_kz);
            
            lambdax_Ts_x = Lambdas.Lambda_x.*Ts_x;
            lambday_Ts_x = Lambdas.Lambda_y.*Ts_x;
            lambdaz_Ts_x = Lambdas.Lambda_z.*Ts_x;
            
            y_1 = FAME_Matrix_Vector_Production_IFFT_Single_General(lambdax_Ts_x,grid_num(1), grid_num(2), grid_num(3), grid_num(1)*grid_num(2)*grid_num(3),Lambdas.D_kx,Lambdas.D_ky,Lambdas.D_kz);
            y_2 = FAME_Matrix_Vector_Production_IFFT_Single_General(lambday_Ts_x,grid_num(1), grid_num(2), grid_num(3), grid_num(1)*grid_num(2)*grid_num(3),Lambdas.D_kx,Lambdas.D_ky,Lambdas.D_kz);
            y_3 = FAME_Matrix_Vector_Production_IFFT_Single_General(lambdaz_Ts_x,grid_num(1), grid_num(2), grid_num(3), grid_num(1)*grid_num(2)*grid_num(3),Lambdas.D_kx,Lambdas.D_ky,Lambdas.D_kz);

            test_x = norm(C_1*x - y_1);
            test_y = norm(C_2*x - y_2);
            test_z = norm(C_3*x - y_3);
            
            x             = rand(3*grid_num(1)*grid_num(2)*grid_num(3),1);
            Qrs_x         = FAME_Matrix_Vector_Production_Qr_General(x, grid_num(1), grid_num(2), grid_num(3), grid_num(1)*grid_num(2)*grid_num(3), Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_kx, Lambdas.D_ky, Lambdas.D_kz, Lambdas.F_ky, Lambdas.F_kz, 'hermitian');
            Sigma_r_Qrs_x = (Lambdas.Sigma_r.^2).*Qrs_x;
            y             = FAME_Matrix_Vector_Production_Qr_General(Sigma_r_Qrs_x, grid_num(1), grid_num(2), grid_num(3), grid_num(1)*grid_num(2)*grid_num(3), Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_kx, Lambdas.D_ky, Lambdas.D_kz, Lambdas.F_ky, Lambdas.F_kz, 'normal');
    end

    test_CsC = norm(Cs*(C*x) - y);
    if test_x > 1e-2
        flag_C1  = 0;
        flag_all = 0;
    end
    if test_y > 1e-2
        flag_C2  = 0;
        flag_all = 0;
    end
    if test_z > 1e-2
        flag_C3  = 0;
        flag_all = 0;
    end
    if test_CsC > 1e-2
        flag_CsC = 0;
        flag_all = 0;
    end
    flag = [flag_all, flag_C1, flag_C2, flag_C3, flag_CsC ];
    if flag_all ~= 1
        for i = 2:length(flag)
            if flag(i) == 0
                error('The eigendecomposition of C%d occure some bugs. Please let us know the example of this error, we will extremely grateful. \nOur email address is "jiawei.am05g@g2.nctu.edu.tw". Thank you very much. \n', i-1);
            end
        end
    end

end