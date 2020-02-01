function FAME_Tools_CPU2GPU_invB()
    clear
    clc
    
    %% Set parameters
    FAME_Folder_Manager( )

    [ Popt ] = FAME_User_Option();

    grid_nums = [Popt.mesh.x_grid_num Popt.mesh.y_grid_num Popt.mesh.z_grid_num];

    [ Popt.lattice.lattice_vec_a, Popt.mesh.edge_lens, Popt.mesh.mesh_lens, Omega, Popt.lattice.lattice_type, Popt.lattice.lattice_constant, Popt.material ] = ...
        FAME_Parameter_Generator( grid_nums, Popt.lattice.lattice_type, Popt.lattice.lattice_constant, Popt.material );
    
    [ ~, wave_vec_num, path_string, ~ ] = ...
        FAME_Parameter_Brillouin_Zone_Path( Popt.lattice.part_num, Popt.lattice.lattice_type, Popt.lattice.lattice_constant, Popt.lattice.lattice_vec_a );
    
    FAME_Tools_Display_Info(grid_nums,Popt.lattice,path_string,Popt.material,wave_vec_num);
    
    %% Set plot axes
    plot_flag = 0;
    if plot_flag
        fig_prim = figure(1);clf
        set(fig_prim,'name','Material shape in primitive cell');
        hax_prim = axes(fig_prim);
        fig_unit = figure(2);clf
        set(fig_unit,'name','Material shape in unit cell');
        hax_unit = axes(fig_unit);

        FAME_Plot_Material_User_Defined(Popt.lattice.lattice_vec_a(:,1),Popt.lattice.lattice_vec_a(:,2),Popt.lattice.lattice_vec_a(:,3), lattice_type, lattice_constant, material.parameters, hax_prim, 'primitive_cell');
        FAME_Plot_Material_User_Defined(Popt.lattice.lattice_vec_a(:,1),Popt.lattice.lattice_vec_a(:,2),Popt.lattice.lattice_vec_a(:,3), lattice_type, lattice_constant, material.parameters, hax_unit, 'unit_cell');
    end
    %% Generate B matrix
    display_material = 'on';
    fprintf('Construct material dependent matrix B......')
    switch Popt.material.material_type
        case 'isotropic'
        [ B_eps, B_mu, invB_eps, invB_mu ] = ...
            FAME_Matrix_B_Isotropic( Popt.material.ele_permitt_out, Popt.material.ele_permitt, Popt.material.mag_permeab, grid_nums, Popt.mesh.mesh_lens, Popt.lattice.lattice_vec_a, Popt.lattice.lattice_constant, Omega, Popt.material.material_handle, display_material, []);
%             ele_permitt_out, ele_permitt_in, mag_permeab, grid_nums, mesh_lens, lattice_vec_a, lattice_constant, Omega, material_handle, display_material, color_map
        case 'anisotropic'
        [ B_eps, B_mu, invB_eps, invB_mu ] = ...
            FAME_Matrix_B_Anisotropic( Popt.material.ele_permitt, Popt.material.mag_permeab, grid_nums, Popt.mesh.mesh_lens, Popt.lattice.lattice_type, Popt.lattice.lattice_vec_a, Popt.lattice.lattice_constant, Omega, Popt.material.material_handle, display_material, [] );
        case 'biisotropic'
        [  B_eps,  B_mu,  B_xi,  B_zeta  ] = ...
            FAME_Matrix_B_Biisotropic( Popt.material.ele_permitt, Popt.material.mag_permeab, Popt.material.reciprocity, Popt.material.chirality, grid_nums, Popt.mesh.mesh_lens, Popt.lattice.lattice_type, Popt.lattice.lattice_vec_a, Popt.lattice.lattice_constant, Omega, Popt.material.material_handle, display_material, [] );    
    end
    % Save as invB.txt file
    fileID = fopen('invB_Hexagonal_n128_p12.txt','w');
    fprintf(fileID,'%f \n',invB_eps);
    fclose(fileID); 

    %% Generate parameter file for C code
    n = grid_nums;
    part_num = Popt.lattice.part_num;
    if strcmp(Popt.lattice.lattice_type,'body_centered_cubic')==1 
        Popt.lattice.lattice_type='BCC';
    end
    if strcmp(Popt.lattice.lattice_type,'face_centered_cubic')==1 
        Popt.lattice.lattice_type='FCC';
    end
    lattice = Popt.lattice.lattice_type;
    path = path_string(1:end-1);
    a     = Popt.lattice.lattice_constant.a;
    b     = Popt.lattice.lattice_constant.b;
    c     = Popt.lattice.lattice_constant.c;
    alpha = Popt.lattice.lattice_constant.alpha;
    beta  = Popt.lattice.lattice_constant.beta;
    gamma = Popt.lattice.lattice_constant.gamma;

    fid=fopen('PC_para.txt','wt');
    str = '# aaa';
    fprintf(fid,'%s \n',str);
    AA = num2str(n);
    fprintf(fid,'%s \n',AA);
    fprintf(fid,'%s \n',str);
    AA = num2str(part_num);
    fprintf(fid,'%s \n',AA);
    fprintf(fid,'%s \n',str);
    AA = '1.0';
    fprintf(fid,'%s \n',AA);
    fprintf(fid,'%s \n',str);
    AA = lattice;
    fprintf(fid,'%s \n',AA);
    fprintf(fid,'%s \n',str);
    AA = 'iso Use';
    fprintf(fid,'%s \n',AA);
    fprintf(fid,'%s \n',str);
    AA = path;
    fprintf(fid,'%s \n',AA);
    fprintf(fid,'%s \n',str);
    AA = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18';
    fprintf(fid,'%s \n',AA);
    fprintf(fid,'%s \n',str);
    AA = '1.0 0.0';
    fprintf(fid,'%s \n',AA);
    fprintf(fid,'%s \n',str);
    fprintf(fid,'%s \n',str);
    AA = '0.12';
    fprintf(fid,'%s \n',AA);
    fprintf(fid,'%s \n',str);
    AA = '0.11';
    fprintf(fid,'%s \n',AA);
    fprintf(fid,'%s \n',str);
    AA = 'User';
    fprintf(fid,'%s \n',AA);
    fprintf(fid,'%s \n',str);
    fprintf(fid,'%s \n',str);
    AA = '1.1';
    fprintf(fid,'%s \n',AA);
    fprintf(fid,'%s \n',str);
    AA = '10';
    fprintf(fid,'%s \n',AA);
    fprintf(fid,'%s \n',str);
    AA = '0.0';
    fprintf(fid,'%s \n',AA);
    fprintf(fid,'%s \n',str);
    AA = '0.25 -0.125 0.5';
    fprintf(fid,'%s \n',AA);
    fprintf(fid,'%s \n',str);
    AA = num2str([a b c alpha beta gamma]);
    fprintf(fid,'%s \n',AA);
    fclose(fid);
end