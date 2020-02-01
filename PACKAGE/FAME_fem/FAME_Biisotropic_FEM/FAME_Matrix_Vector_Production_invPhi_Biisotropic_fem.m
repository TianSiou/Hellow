function vec_y_ele = FAME_Matrix_Vector_Production_invPhi_Biisotropic_fem(temp_vec_y_ele, B)
    temp_vec_y_ele_1 = temp_vec_y_ele(1:end/3);
    temp_vec_y_ele_2 = temp_vec_y_ele(end/3+1:2*end/3);
    temp_vec_y_ele_3 = temp_vec_y_ele(2*end/3+1:end);
    if strcmp(B.lssvr, 'mldivide')
        vec_y_ele = B.Phi\temp_vec_y_ele;
    elseif strcmp(B.lssvr, 'pcg')
        [ vec_y_ele, ~, ~, ~  ] = pcg( B.Phi, temp_vec_y_ele, 1e-12, 1000);
%         [ vec_y_ele, ~, ~, inner_inner_iter(inner_inner_count) ] = pcg( B.Phi, temp_vec_y_ele, 1e-12, 1000);
%         inner_inner_cpu_time(inner_inner_count) = toc(time_start_pcg_B);
    elseif isfield(B,'P_1')
        temp_vec_y_ele_1 = B.L_1\temp_vec_y_ele_1(B.P_1);
        temp_vec_y_ele_1 = B.U_1\temp_vec_y_ele_1;
        vec_y_ele_1      = temp_vec_y_ele_1(B.invP_1);
        temp_vec_y_ele_2 = B.L_2\temp_vec_y_ele_2(B.P_2);
        temp_vec_y_ele_2 = B.U_2\temp_vec_y_ele_2;
        vec_y_ele_2      = temp_vec_y_ele_2(B.invP_2);
        temp_vec_y_ele_3 = B.L_3\temp_vec_y_ele_3(B.P_3);
        temp_vec_y_ele_3 = B.U_3\temp_vec_y_ele_3;
        vec_y_ele_3      = temp_vec_y_ele_3(B.invP_3);
        
        vec_y_ele = [vec_y_ele_1;vec_y_ele_2;vec_y_ele_3];
    else
        temp_vec_y_ele_1 = B.L_1\temp_vec_y_ele_1;
        vec_y_ele_1      = B.U_1\temp_vec_y_ele_1;
        temp_vec_y_ele_2 = B.L_2\temp_vec_y_ele_2;
        vec_y_ele_2      = B.U_2\temp_vec_y_ele_2;
        temp_vec_y_ele_3 = B.L_3\temp_vec_y_ele_3;
        vec_y_ele_3      = B.U_3\temp_vec_y_ele_3;
        
        vec_y_ele = [vec_y_ele_1;vec_y_ele_2;vec_y_ele_3];
    end
end