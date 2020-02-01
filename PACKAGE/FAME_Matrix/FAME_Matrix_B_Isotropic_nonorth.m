function B = FAME_Matrix_B_Isotropic_nonorth( Par_mesh, Par_material)

N = Par_mesh.grid_num(1)*Par_mesh.grid_num(2)*Par_mesh.grid_num(3);
N_material_handle = length(Par_material.B.org_idx);
N_permitt         = length(Par_material.ele_permitt_in);
% N_permeab         = length(Par_material.mag_permeab_in);

if (N_material_handle~=N_permitt)
    fprintf('\n');
    warning('The input number of material handle, permittivity and permeability not equal! Please check these input data.')
    ele_permitt_in = Par_material.ele_permitt_in(1)*ones(N_material_handle,1);
else
    ele_permitt_in = Par_material.ele_permitt_in;
end

B_eps_x = Par_material.ele_permitt_out*ones(N,1);  B_eps_y = Par_material.ele_permitt_out*ones(N,1);  B_eps_z = Par_material.ele_permitt_out*ones(N,1);
for i = 1:N_material_handle
    B_eps_x(Par_material.B.org_idx{i}) = ele_permitt_in(i);
    B_eps_y(Par_material.B.org_idx{i}) = ele_permitt_in(i);
    B_eps_z(Par_material.B.org_idx{i}) = ele_permitt_in(i);
end
B.B_eps    = [ B_eps_x; B_eps_y; B_eps_z ];
B.invB_eps = 1./B.B_eps;
end
