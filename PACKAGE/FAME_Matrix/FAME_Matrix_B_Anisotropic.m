function B = FAME_Matrix_B_Anisotropic( Par_mesh, Par_material )

N = Par_mesh.grid_num(1)*Par_mesh.grid_num(2)*Par_mesh.grid_num(3);
N_material_handle = length(Par_material.B.ele_x_idx);
N_permitt         = length(Par_material.ele_permitt_in);
N_permeab         = length(Par_material.mag_permeab_in);

if (N_material_handle~=N_permitt) || (N_material_handle~=N_permeab) || (N_permitt~=N_permeab)
    error('The input number of material handle, permittivity and permeability not equal! Please check these input data.')
end

B_eps_x = ones(N,1);  B_eps_y = ones(N,1);  B_eps_z = ones(N,1);
B_mu_x  = ones(N,1);  B_mu_y  = ones(N,1);  B_mu_z  = ones(N,1);
for i = 1:N_material_handle
%     ele_permitt_i     = Omega*Par_material.ele_permitt_in{i}*Omega';
%     mag_permeab_i     = Omega*Par_material.mag_permeab_in{i}*Omega';
    ele_permitt_i     = Par_material.ele_permitt_in{i};
    mag_permeab_i     = Par_material.mag_permeab_in{i};

    [ ele_permitt_rotate{i}, ele_permitt_diag ] = eig(ele_permitt_i);
    [ mag_permeab_rotate{i}, mag_permeab_diag ] = eig(mag_permeab_i);

    B_eps_x(Par_material.B.ele_x_idx{i}) = ele_permitt_diag(1,1);
    B_eps_y(Par_material.B.ele_y_idx{i}) = ele_permitt_diag(2,2);
    B_eps_z(Par_material.B.ele_z_idx{i}) = ele_permitt_diag(3,3);
    B_mu_x(Par_material.B.mag_x_idx{i})  = mag_permeab_diag(1,1);
    B_mu_y(Par_material.B.mag_y_idx{i})  = mag_permeab_diag(2,2);
    B_mu_z(Par_material.B.mag_z_idx{i})  = mag_permeab_diag(3,3);
    
%     I    = speye(N);
%     Q_ele_permitt    = [ ele_permitt_rotate(1,1)*I,ele_permitt_rotate(1,2)*I,ele_permitt_rotate(1,3)*I;  
%                          ele_permitt_rotate(2,1)*I,ele_permitt_rotate(2,2)*I,ele_permitt_rotate(2,3)*I;
%                          ele_permitt_rotate(3,1)*I,ele_permitt_rotate(3,2)*I,ele_permitt_rotate(3,3)*I];
%     Q_mag_permeab    = [ mag_permeab_rotate(1,1)*I,mag_permeab_rotate(1,2)*I,mag_permeab_rotate(1,3)*I;  
%                          mag_permeab_rotate(2,1)*I,mag_permeab_rotate(2,2)*I,mag_permeab_rotate(2,3)*I;
%                          mag_permeab_rotate(3,1)*I,mag_permeab_rotate(3,2)*I,mag_permeab_rotate(3,3)*I];                 
end
B.B_eps    = blk_33_inverse(B_eps_x, B_eps_y, B_eps_z, ele_permitt_rotate, Par_material.B.ele_x_idx, N);
B.B_mu     = blk_33_inverse(B_mu_x , B_mu_y , B_mu_z , mag_permeab_rotate, Par_material.B.mag_x_idx, N);
B.invB_eps = blk_33_inverse(1./B_eps_x, 1./B_eps_y, 1./B_eps_z, ele_permitt_rotate, Par_material.B.ele_x_idx, N);
B.invB_mu  = blk_33_inverse(1./B_mu_x ,  1./B_mu_y, 1./B_mu_z , mag_permeab_rotate, Par_material.B.mag_x_idx, N);

% B_eps    = Q_ele_permitt * spdiags( [B_eps_x; B_eps_y; B_eps_z]   , 0, 3*N, 3*N ) * Q_ele_permitt';
% B_mu     = Q_mag_permeab * spdiags( [B_mu_x; B_mu_y; B_mu_z]      , 0, 3*N, 3*N ) * Q_mag_permeab';
% invB_eps = Q_ele_permitt * spdiags( 1./[B_eps_x; B_eps_y; B_eps_z], 0, 3*N, 3*N ) * Q_ele_permitt';
% invB_mu  = Q_mag_permeab * spdiags( 1./[B_mu_x; B_mu_y; B_mu_z]   , 0, 3*N, 3*N ) * Q_mag_permeab';
end

function B = blk_33_inverse(Bx,By,Bz,Q_cell,idx,N)
B.d11 = zeros(N,1); B.d12 = zeros(N,1); B.d13 = zeros(N,1);
B.d21 = zeros(N,1); B.d22 = zeros(N,1); B.d23 = zeros(N,1);
B.d31 = zeros(N,1); B.d32 = zeros(N,1); B.d33 = zeros(N,1);
for i = 1:length(idx)
    temp_idx = idx{i};
    Q        = Q_cell{i};

    B.d11(temp_idx) = Q(1,1)*conj(Q(1,1)) * Bx(temp_idx) + Q(1,2)*conj(Q(1,2)) * By(temp_idx) + Q(1,3)*conj(Q(1,3)) * Bz(temp_idx);
    B.d12(temp_idx) = Q(1,1)*conj(Q(2,1)) * Bx(temp_idx) + Q(1,2)*conj(Q(2,2)) * By(temp_idx) + Q(1,3)*conj(Q(2,3)) * Bz(temp_idx);
    B.d13(temp_idx) = Q(1,1)*conj(Q(3,1)) * Bx(temp_idx) + Q(1,2)*conj(Q(3,2)) * By(temp_idx) + Q(1,3)*conj(Q(3,3)) * Bz(temp_idx);
    B.d21(temp_idx) = Q(2,1)*conj(Q(1,1)) * Bx(temp_idx) + Q(2,2)*conj(Q(1,2)) * By(temp_idx) + Q(2,3)*conj(Q(1,3)) * Bz(temp_idx);
    B.d22(temp_idx) = Q(2,1)*conj(Q(2,1)) * Bx(temp_idx) + Q(2,2)*conj(Q(2,2)) * By(temp_idx) + Q(2,3)*conj(Q(2,3)) * Bz(temp_idx);
    B.d23(temp_idx) = Q(2,1)*conj(Q(3,1)) * Bx(temp_idx) + Q(2,2)*conj(Q(3,2)) * By(temp_idx) + Q(2,3)*conj(Q(3,3)) * Bz(temp_idx);
    B.d31(temp_idx) = Q(3,1)*conj(Q(1,1)) * Bx(temp_idx) + Q(3,2)*conj(Q(1,2)) * By(temp_idx) + Q(3,3)*conj(Q(1,3)) * Bz(temp_idx);
    B.d32(temp_idx) = Q(3,1)*conj(Q(2,1)) * Bx(temp_idx) + Q(3,2)*conj(Q(2,2)) * By(temp_idx) + Q(3,3)*conj(Q(2,3)) * Bz(temp_idx);
    B.d33(temp_idx) = Q(3,1)*conj(Q(3,1)) * Bx(temp_idx) + Q(3,2)*conj(Q(3,2)) * By(temp_idx) + Q(3,3)*conj(Q(3,3)) * Bz(temp_idx);
    rest_idx = setdiff(1:N,temp_idx);
end
B.d11(rest_idx) = Q(1,1)*conj(Q(1,1)) * Bx(rest_idx) + Q(1,2)*conj(Q(1,2)) * By(rest_idx) + Q(1,3)*conj(Q(1,3)) * Bz(rest_idx);
B.d12(rest_idx) = Q(1,1)*conj(Q(2,1)) * Bx(rest_idx) + Q(1,2)*conj(Q(2,2)) * By(rest_idx) + Q(1,3)*conj(Q(2,3)) * Bz(rest_idx);
B.d13(rest_idx) = Q(1,1)*conj(Q(3,1)) * Bx(rest_idx) + Q(1,2)*conj(Q(3,2)) * By(rest_idx) + Q(1,3)*conj(Q(3,3)) * Bz(rest_idx);
B.d21(rest_idx) = Q(2,1)*conj(Q(1,1)) * Bx(rest_idx) + Q(2,2)*conj(Q(1,2)) * By(rest_idx) + Q(2,3)*conj(Q(1,3)) * Bz(rest_idx);
B.d22(rest_idx) = Q(2,1)*conj(Q(2,1)) * Bx(rest_idx) + Q(2,2)*conj(Q(2,2)) * By(rest_idx) + Q(2,3)*conj(Q(2,3)) * Bz(rest_idx);
B.d23(rest_idx) = Q(2,1)*conj(Q(3,1)) * Bx(rest_idx) + Q(2,2)*conj(Q(3,2)) * By(rest_idx) + Q(2,3)*conj(Q(3,3)) * Bz(rest_idx);
B.d31(rest_idx) = Q(3,1)*conj(Q(1,1)) * Bx(rest_idx) + Q(3,2)*conj(Q(1,2)) * By(rest_idx) + Q(3,3)*conj(Q(1,3)) * Bz(rest_idx);
B.d32(rest_idx) = Q(3,1)*conj(Q(2,1)) * Bx(rest_idx) + Q(3,2)*conj(Q(2,2)) * By(rest_idx) + Q(3,3)*conj(Q(2,3)) * Bz(rest_idx);
B.d33(rest_idx) = Q(3,1)*conj(Q(3,1)) * Bx(rest_idx) + Q(3,2)*conj(Q(3,2)) * By(rest_idx) + Q(3,3)*conj(Q(3,3)) * Bz(rest_idx);
end
