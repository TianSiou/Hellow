clear all
clc

FAME_Folder_Manager([]) % loading path for whole FAME

%% Load user options from exist file "FAME_User_Option.m"
% [ Popt ] = FAME_User_Option();
[ Popt ] = FAME_User_Option();
%% Generate modified lattice vectors and lattice constants for computing
[ Par.mesh, Par.lattice, Par.recip_lattice, Par.material, Par.eig ] = FAME_Parameter_Generator( Popt );
%% Generate wave vector array by partition the path of Brillouin zone 
[ Par.recip_lattice ] = FAME_Parameter_Brillouin_Zone_Path( Popt.recip_lattice.part_num, Par.lattice, Par.recip_lattice );
%% Locating indices for the material inside
[  Par.material.B.org_idx ] = FAME_Material_Locate_Index_nonorth( Par.mesh, Par.lattice, Par.material);
%% Start FAME
FAME_option.discrete_method = 'Yee_scheme'; % 'fem' or 'Yee_scheme'
[ omega_array ] = FAME_Main_Code_nonorth( Par.mesh, Par.lattice, Par.recip_lattice, Par.material, Par.eig, Par.recip_lattice.wave_vec_array, FAME_option );
%% Plot band structure
FAME_Plot_Band_Structure( Par.recip_lattice.path_string_new, Par.recip_lattice.part_num, omega_array/(2*pi) )