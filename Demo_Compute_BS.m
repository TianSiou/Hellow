  clear all
clc
FAME_Folder_Manager([]) % loading path for whole FAME

% mean_iter = [];
% restart_cout = [];
% cpu_time = [];

%  for i = 1 : 1
%      for j = 1 : 5
         clear comput_info ele_field FAME_option mag_field n_supercell omega_array Par Popt 

%% Load user options from exist file "FAME_User_Option.m"
[ Popt ] = FAME_User_Option();

% n_supercell = [j j];
% Popt.mesh.grid_num = [16 * i, 16 * i, 16 * i].*[sum(n_supercell), 1, 1]; % The grid numbers
% Popt.lattice.lattice_vec_a = [ 0.5*[-1;1;1], sum(n_supercell)*0.5*[1;-1;1], 0.5*[1;1;-1]];
%% Generate modified lattice vectors and lattice constants for computing
[ Par.mesh, Par.lattice, Par.recip_lattice, Par.material, Par.eig ] = FAME_Parameter_Generator( Popt );
%% Generate wave vector array by partition the path of Brillouin zone 
path_str = 'G';
% path_str = 'GXMGRX|MR';
[ Par.recip_lattice ] = FAME_Parameter_Brillouin_Zone_Path( Popt.recip_lattice.part_num, Par.lattice, Par.recip_lattice, path_str);
%% Locating indices for the material inside
[ Par.material.B.ele_x_idx, Par.material.B.ele_y_idx, Par.material.B.ele_z_idx, Par.material.B.mag_x_idx, Par.material.B.mag_y_idx, Par.material.B.mag_z_idx, Par.material.B.org_idx ] = ...
    FAME_Material_Locate_Index( Par.mesh, Par.lattice, Par.material);
%% Start FAME
fprintf('Begin time : %s',datestr(now))
% count = (i-1) * 10 + j
FAME_option.discrete_method = 'Yee_scheme'; % 'fem' or 'Yee_scheme'
[ omega_array, ele_field, mag_field, comput_info ] = FAME_Main_Code( Par.mesh, Par.lattice, Par.material, Par.eig, Par.recip_lattice.wave_vec_array, FAME_option );
%% Plot band structure
% FAME_Plot_Band_Structure( Par.recip_lattice.path_string_new, Par.recip_lattice.part_num, omega_array/(2*pi) )
% eval(['save 1Cyli_1_sup_shift45_16',' -v7.3'])
% eval(['save SDGsup4_1_shift45_Precon_32',' -v7.3'])
% mean_iter(i , j) =  mean(comput_info.LS_iter{1, 1});
% restart_cout(i , j) = length(comput_info.LS_iter{1, 1});
% cpu_time(i , j) = comput_info.cpu_time;
%      end
%  end