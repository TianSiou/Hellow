function B = FAME_Matrix_B_DrudeModel( Par_mesh, Par_material, freq)

    N = Par_mesh.grid_num(1)*Par_mesh.grid_num(2)*Par_mesh.grid_num(3);
    N_material_handle = length(Par_material.B.ele_x_idx);

    B_eps_x = Par_material.ele_permitt_out*ones(N,1);  B_eps_y = Par_material.ele_permitt_out*ones(N,1);  B_eps_z = Par_material.ele_permitt_out*ones(N,1);
    
    ele_permitt_in = Par_material.ele_permitt_in(freq);
    fprintf('(£` = %2.2f + %2.2fi)',real(ele_permitt_in),imag(ele_permitt_in));
    for i = 1:N_material_handle
        B_eps_x(Par_material.B.ele_x_idx{i}) = ele_permitt_in;
        B_eps_y(Par_material.B.ele_y_idx{i}) = ele_permitt_in;
        B_eps_z(Par_material.B.ele_z_idx{i}) = ele_permitt_in;
    end
    B.B_eps    = [ B_eps_x; B_eps_y; B_eps_z ];
    B.invB_eps = 1./B.B_eps;
end

