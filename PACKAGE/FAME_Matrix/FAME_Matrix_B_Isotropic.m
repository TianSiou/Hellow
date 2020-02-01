function B = FAME_Matrix_B_Isotropic( Par_mesh, Par_material)

N = Par_mesh.grid_num(1)*Par_mesh.grid_num(2)*Par_mesh.grid_num(3);
N_material_handle = length(Par_material.B.ele_x_idx);
N_permitt         = length(Par_material.ele_permitt_in);
% N_permeab         = length(Par_material.mag_permeab_in);

if (N_material_handle~=N_permitt)
    fprintf('\n');
    warning('The input number of material handle, permittivity and permeability not equal! Please check these input data.')
    ele_permitt_in = Par_material.ele_permitt_in(1)*ones(N_material_handle,1);
%     mag_permeab_in = Par_material.mag_permeab_in(1)*ones(N_material_handle,1);
else
    ele_permitt_in = Par_material.ele_permitt_in;
%     mag_permeab_in = Par_material.mag_permeab_in;
end

B_eps_x = Par_material.ele_permitt_out*ones(N,1);  B_eps_y = Par_material.ele_permitt_out*ones(N,1);  B_eps_z = Par_material.ele_permitt_out*ones(N,1);
% B_mu_x  = Par_material.mag_permeab_out*ones(N,1);  B_mu_y  = Par_material.mag_permeab_out*ones(N,1);  B_mu_z  = Par_material.mag_permeab_out*ones(N,1);

for i = 1:N_material_handle
    B_inout_ele_1 = zeros(N,1); B_inout_ele_2 = zeros(N,1); B_inout_ele_3 = zeros(N,1);
%     B_inout_mag_1 = zeros(N,1); B_inout_mag_2 = zeros(N,1); B_inout_mag_3 = zeros(N,1);
    B_inout_ele_1(Par_material.B.ele_x_idx{1}) = 1; B_inout_ele_2(Par_material.B.ele_y_idx{1}) = 1; B_inout_ele_3(Par_material.B.ele_z_idx{1}) = 1;
%     B_inout_mag_1(Par_material.B.mag_x_idx{1}) = 1; B_inout_mag_2(Par_material.B.mag_y_idx{1}) = 1; B_inout_mag_3(Par_material.B.mag_z_idx{1}) = 1;
    B.B_inout_ele(:,i) = [B_inout_ele_1;B_inout_ele_2;B_inout_ele_3];
%     B.B_inout_mag(:,i) = [B_inout_mag_1;B_inout_mag_2;B_inout_mag_3];
    
    B_eps_x(Par_material.B.ele_x_idx{i}) = ele_permitt_in(i);
    B_eps_y(Par_material.B.ele_y_idx{i}) = ele_permitt_in(i);
    B_eps_z(Par_material.B.ele_z_idx{i}) = ele_permitt_in(i);
%     B_mu_x(Par_material.B.mag_x_idx{i})  = mag_permeab_in(i);
%     B_mu_y(Par_material.B.mag_y_idx{i})  = mag_permeab_in(i);
%     B_mu_z(Par_material.B.mag_z_idx{i})  = mag_permeab_in(i);
end
B.B_eps    = [ B_eps_x; B_eps_y; B_eps_z ];
% B.B_mu     = [ B_mu_x ; B_mu_y ; B_mu_z  ];
B.invB_eps = 1./B.B_eps;
% B.invB_mu  = 1./B.B_mu;
end
