function B = FAME_Matrix_B_Biisotropic( Par_mesh, Par_material )


N                 = Par_mesh.grid_num(1)*Par_mesh.grid_num(2)*Par_mesh.grid_num(3);
N_material_handle = length(Par_material.material_handle);
N_permitt         = length(Par_material.ele_permitt_in);
N_permeab         = length(Par_material.mag_permeab_in);
N_reciprocity     = length(Par_material.reciprocity_in);
N_chirality       = length(Par_material.chirality_in);

if (N_material_handle~=N_permitt) || (N_material_handle~=N_permeab) || (N_material_handle~=N_reciprocity) || (N_material_handle~=N_chirality)
    error('The input number of material handle, permittivity and permeability not equal! Please check these input data.')
end

B_eps_x   = Par_material.ele_permitt_out*ones(N,1);   B_eps_y   = Par_material.ele_permitt_out*ones(N,1);   B_eps_z   = Par_material.ele_permitt_out*ones(N,1);
B_mu_x    = Par_material.mag_permeab_out*ones(N,1);   B_mu_y    = Par_material.mag_permeab_out*ones(N,1);   B_mu_z    = Par_material.mag_permeab_out*ones(N,1);
B_xi_x    = zeros(N,1);  B_xi_y    = zeros(N,1);  B_xi_z    = zeros(N,1);
B_zeta_x  = zeros(N,1);  B_zeta_y  = zeros(N,1);  B_zeta_z  = zeros(N,1);
for i = 1:N_material_handle
    B_eps_x(Par_material.B.ele_x_idx{i}) = Par_material.ele_permitt_in(i);
    B_eps_y(Par_material.B.ele_y_idx{i}) = Par_material.ele_permitt_in(i);
    B_eps_z(Par_material.B.ele_z_idx{i}) = Par_material.ele_permitt_in(i);
    
    B_mu_x(Par_material.B.mag_x_idx{i})  = Par_material.mag_permeab_in(i);
    B_mu_y(Par_material.B.mag_y_idx{i})  = Par_material.mag_permeab_in(i);
    B_mu_z(Par_material.B.mag_z_idx{i})  = Par_material.mag_permeab_in(i);
    %  (xi)£i = £q - i£e and (zeta)£a = £q + i£e
    B_xi_x(Par_material.B.ele_x_idx{i})    = Par_material.reciprocity_in(i) - 1i*Par_material.chirality_in(i);
    B_xi_y(Par_material.B.ele_y_idx{i})    = Par_material.reciprocity_in(i) - 1i*Par_material.chirality_in(i);
    B_xi_z(Par_material.B.ele_z_idx{i})    = Par_material.reciprocity_in(i) - 1i*Par_material.chirality_in(i);
    
    B_zeta_x(Par_material.B.ele_x_idx{i})  = Par_material.reciprocity_in(i) + 1i*Par_material.chirality_in(i);
    B_zeta_y(Par_material.B.ele_y_idx{i})  = Par_material.reciprocity_in(i) + 1i*Par_material.chirality_in(i);
    B_zeta_z(Par_material.B.ele_z_idx{i})  = Par_material.reciprocity_in(i) + 1i*Par_material.chirality_in(i);
end
B.B_eps  = [B_eps_x;B_eps_y;B_eps_z];
B.B_mu   = [B_mu_x;B_mu_y;B_mu_z];
B.B_xi   = [B_xi_x;B_xi_y;B_xi_z];
B.B_zeta = [B_zeta_x;B_zeta_y;B_zeta_z];
end