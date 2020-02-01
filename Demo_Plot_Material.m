clear all
clc

FAME_Folder_Manager([]) % loading path for whole FAME

%% Load user options from exist file "FAME_User_Option.m"
[ Popt ] = FAME_User_Option();
%% Generate modified lattice vectors and lattice constants for computing
[ Par.mesh, Par.lattice, Par.recip_lattice, Par.material, Par.eig ] = FAME_Parameter_Generator( Popt );
%% Locating indices for the material inside
% fig_grid = figure('Grid inner material');clf
% hax_grid = axes(fig_grid);
[ Par.material.B.ele_x_idx, Par.material.B.ele_y_idx, Par.material.B.ele_z_idx, Par.material.B.mag_x_idx, Par.material.B.mag_y_idx, Par.material.B.mag_z_idx, Par.material.B.org_idx ] = ...
    FAME_Material_Locate_Index( Par.mesh, Par.lattice, Par.material);
%% Plot Material shape
figure('Name','Material shape: Primitive cell'); clf
hax_material_primitive = cla;
plot_mode = 'primitive_cell';
FAME_Plot_Material_User_Defined(Par.mesh.grid_num, Par.lattice.lattice_vec_a(:,1), Par.lattice.lattice_vec_a(:,2), Par.lattice.lattice_vec_a(:,3), Par.lattice.lattice_type, Par.lattice.lattice_constant, Par.lattice.Omega, Par.material, hax_material_primitive, plot_mode);
figure('Name','Material shape: Unit cell'); clf
hax_material_unit   = cla;
plot_mode = 'unit_cell';
FAME_Plot_Material_User_Defined(Par.mesh.grid_num, Par.lattice.lattice_vec_a(:,1), Par.lattice.lattice_vec_a(:,2), Par.lattice.lattice_vec_a(:,3), Par.lattice.lattice_type, Par.lattice.lattice_constant, Par.lattice.Omega, Par.material, hax_material_unit, plot_mode);
figure('Name','Material shape: Computational cell'); clf
hax_material_unit   = cla;
plot_mode = 'computational_cell';
FAME_Plot_Material_User_Defined(Par.mesh.grid_num, Par.lattice.lattice_vec_a(:,1), Par.lattice.lattice_vec_a(:,2), Par.lattice.lattice_vec_a(:,3), Par.lattice.lattice_type, Par.lattice.lattice_constant, Par.lattice.Omega, Par.material, hax_material_unit, plot_mode);
% FAME_Plot_Material_User_Defined(grid_nums, Par.lattice.lattice_vec_a(:,1), Par.lattice.lattice_vec_a(:,2), Par.lattice.lattice_vec_a(:,3), Par.lattice.lattice_type, Par.lattice.lattice_constant, Par.lattice.Omega, Par.material, hax_material, plot_mode);