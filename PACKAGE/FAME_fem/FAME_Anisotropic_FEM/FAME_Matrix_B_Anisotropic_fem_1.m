function [ B ] = FAME_Matrix_B_Anisotropic_fem_1( Par_mesh, Par_lattice, Par_material  )

eps_in   = Par_material.ele_permitt_in{1};
eps_out  = Par_material.ele_permitt_out{1};
inner_idx = Par_material.B.org_idx{1};

%% Preprocess
N = Par_mesh.grid_num(1)*Par_mesh.grid_num(2)*Par_mesh.grid_num(3);

Ones_x = ones(Par_mesh.grid_num(1),1);
Ones_y = ones(Par_mesh.grid_num(2),1);
Ones_z = ones(Par_mesh.grid_num(3),1);

Ix  = speye(Par_mesh.grid_num(1));  Iy  = speye(Par_mesh.grid_num(2));  Iz  = speye(Par_mesh.grid_num(3));

B_eps_11 = eps_out(1,1)*ones(N,1);    B_eps_11(inner_idx)  = eps_in(1,1);
B_eps_12 = zeros(N,1);                B_eps_12(inner_idx)  = eps_in(1,2);
B_eps_13 = zeros(N,1);                B_eps_13(inner_idx)  = eps_in(1,3);
B_eps_21 = zeros(N,1);                B_eps_21(inner_idx)  = eps_in(2,1);
B_eps_22 = eps_out(2,2)*ones(N,1);    B_eps_22(inner_idx)  = eps_in(2,2);
B_eps_23 = zeros(N,1);                B_eps_23(inner_idx)  = eps_in(2,3);
B_eps_31 = zeros(N,1);                B_eps_31(inner_idx)  = eps_in(3,1);
B_eps_32 = zeros(N,1);                B_eps_32(inner_idx)  = eps_in(3,2);
B_eps_33 = eps_out(3,3)*ones(N,1);    B_eps_33(inner_idx)  = eps_in(3,3);

%% Construct Matrices
switch Par_lattice.lattice_type
    case 'simple_cubic'
        Kx_forward_noquasi = spdiags(Ones_x,1,Par_mesh.grid_num(1),Par_mesh.grid_num(1));
        Kx_forward_noquasi(Par_mesh.grid_num(1),1) = 1;
        Ky_forward_noquasi = spdiags(Ones_y,1,Par_mesh.grid_num(2),Par_mesh.grid_num(2));
        Ky_forward_noquasi(Par_mesh.grid_num(2),1) = 1;
        Kz_forward_noquasi = spdiags(Ones_z,1,Par_mesh.grid_num(3),Par_mesh.grid_num(3));
        Kz_forward_noquasi(Par_mesh.grid_num(3),1) = 1;

        B_eps_11_shift_x = kron(Iz,kron(Iy,Kx_forward_noquasi))*B_eps_11;
        B_eps_12_shift_x = kron(Iz,kron(Iy,Kx_forward_noquasi))*B_eps_12;
        B_eps_13_shift_x = kron(Iz,kron(Iy,Kx_forward_noquasi))*B_eps_13;

        B_eps_21_shift_y = kron(Iz,kron(Ky_forward_noquasi,Ix))*B_eps_21;
        B_eps_22_shift_y = kron(Iz,kron(Ky_forward_noquasi,Ix))*B_eps_22;
        B_eps_23_shift_y = kron(Iz,kron(Ky_forward_noquasi,Ix))*B_eps_23;

        B_eps_31_shift_z = kron(Kz_forward_noquasi,kron(Iy,Ix))*B_eps_31;
        B_eps_32_shift_z = kron(Kz_forward_noquasi,kron(Iy,Ix))*B_eps_32;
        B_eps_33_shift_z = kron(Kz_forward_noquasi,kron(Iy,Ix))*B_eps_33;
    case 'face_centered_cubic'
        J2_noquasi =   spdiags([ Ones_x, exp(0)*Ones_x ], [-0.5*Par_mesh.grid_num(1),0.5*Par_mesh.grid_num(1)], Par_mesh.grid_num(1),Par_mesh.grid_num(1));
        J3_noquasi = [ sparse((Par_mesh.grid_num(2)*Par_mesh.grid_num(1))/3,2*(Par_mesh.grid_num(2)*Par_mesh.grid_num(1)/3)),    exp(0)*speye((Par_mesh.grid_num(2)*Par_mesh.grid_num(1)/3)) ;...
                       kron(speye(2*Par_mesh.grid_num(2)/3),J2_noquasi)                                     ,    sparse(2*(Par_mesh.grid_num(2)*Par_mesh.grid_num(1)/3),(Par_mesh.grid_num(2)*Par_mesh.grid_num(1))/3) ];
        
        Kx_forward_noquasi = spdiags(Ones_x,1,Par_mesh.grid_num(1),Par_mesh.grid_num(1));
        Kx_forward_noquasi(Par_mesh.grid_num(1),1) = exp(0);
        Kx_forward_noquasi = kron(Iz,kron(Iy,Kx_forward_noquasi));
        Ky_forward_noquasi = spdiags(Ones_y,1,Par_mesh.grid_num(2),Par_mesh.grid_num(2));
        Ky_forward_noquasi = kron(Ky_forward_noquasi,Ix);
        Ky_forward_noquasi((Par_mesh.grid_num(2)-1)*Par_mesh.grid_num(1)+1:end, 1:Par_mesh.grid_num(1)) = J2_noquasi;
        Ky_forward_noquasi = kron(Iz,Ky_forward_noquasi);
        Kz_forward_noquasi = spdiags(Ones_z,1,Par_mesh.grid_num(3),Par_mesh.grid_num(3));
        Kz_forward_noquasi = kron(Kz_forward_noquasi,kron(Iy,Ix));
        Kz_forward_noquasi((Par_mesh.grid_num(3)-1)*Par_mesh.grid_num(2)*Par_mesh.grid_num(1)+1:end, 1:Par_mesh.grid_num(2)*Par_mesh.grid_num(1)) = J3_noquasi;
                   
        B_eps_11_shift_x = Kx_forward_noquasi*B_eps_11;
        B_eps_12_shift_x = Kx_forward_noquasi*B_eps_12;
        B_eps_13_shift_x = Kx_forward_noquasi*B_eps_13;

        B_eps_21_shift_y = Ky_forward_noquasi*B_eps_21;
        B_eps_22_shift_y = Ky_forward_noquasi*B_eps_22;
        B_eps_23_shift_y = Ky_forward_noquasi*B_eps_23;

        B_eps_31_shift_z = Kz_forward_noquasi*B_eps_31;
        B_eps_32_shift_z = Kz_forward_noquasi*B_eps_32;
        B_eps_33_shift_z = Kz_forward_noquasi*B_eps_33;
    case 'body_centered_cubic'
        J1_noquasi = spdiags([ Ones_x, Ones_x ], [1/3*Par_mesh.grid_num(1),-2/3*Par_mesh.grid_num(1)], Par_mesh.grid_num(1),Par_mesh.grid_num(1));
        J2_noquasi = spdiags([ Ones_x, Ones_x ], [2/3*Par_mesh.grid_num(1),-1/3*Par_mesh.grid_num(1)], Par_mesh.grid_num(1),Par_mesh.grid_num(1));
        J3_noquasi = [ sparse(1/2*Par_mesh.grid_num(2)*Par_mesh.grid_num(1),1/2*Par_mesh.grid_num(2)*Par_mesh.grid_num(1)) , sparse(kron(speye(0.5*Par_mesh.grid_num(2)), J1_noquasi)) ;...
                       sparse(kron(speye(0.5*Par_mesh.grid_num(2)),J2_noquasi)), sparse(1/2*Par_mesh.grid_num(2)*Par_mesh.grid_num(1),1/2*Par_mesh.grid_num(2)*Par_mesh.grid_num(1)) ];
        
        Kx_forward_noquasi = spdiags(Ones_x,1,Par_mesh.grid_num(1),Par_mesh.grid_num(1));
        Kx_forward_noquasi(Par_mesh.grid_num(1),1) = exp(0);
        Kx_forward_noquasi = kron(Iz,kron(Iy,Kx_forward_noquasi));
        Ky_forward_noquasi = spdiags(Ones_y,1,Par_mesh.grid_num(2),Par_mesh.grid_num(2));
        Ky_forward_noquasi = kron(Ky_forward_noquasi,Ix);
        Ky_forward_noquasi((Par_mesh.grid_num(2)-1)*Par_mesh.grid_num(1)+1:end, 1:Par_mesh.grid_num(1)) = J1_noquasi;
        Ky_forward_noquasi = kron(Iz,Ky_forward_noquasi);
        Kz_forward_noquasi = spdiags(Ones_z,1,Par_mesh.grid_num(3),Par_mesh.grid_num(3));
        Kz_forward_noquasi = kron(Kz_forward_noquasi,kron(Iy,Ix));
        Kz_forward_noquasi((Par_mesh.grid_num(3)-1)*Par_mesh.grid_num(2)*Par_mesh.grid_num(1)+1:end, 1:Par_mesh.grid_num(2)*Par_mesh.grid_num(1)) = J3_noquasi;
                   
        B_eps_11_shift_x = Kx_forward_noquasi*B_eps_11;
        B_eps_12_shift_x = Kx_forward_noquasi*B_eps_12;
        B_eps_13_shift_x = Kx_forward_noquasi*B_eps_13;

        B_eps_21_shift_y = Ky_forward_noquasi*B_eps_21;
        B_eps_22_shift_y = Ky_forward_noquasi*B_eps_22;
        B_eps_23_shift_y = Ky_forward_noquasi*B_eps_23;

        B_eps_31_shift_z = Kz_forward_noquasi*B_eps_31;
        B_eps_32_shift_z = Kz_forward_noquasi*B_eps_32;
        B_eps_33_shift_z = Kz_forward_noquasi*B_eps_33;
    otherwise
        J1_noquasi = spdiags([ Ones_x, Ones_x ], [1/3*Par_mesh.grid_num(1),-2/3*Par_mesh.grid_num(1)], Par_mesh.grid_num(1),Par_mesh.grid_num(1));
        J2_noquasi = spdiags([ Ones_x, Ones_x ], [2/3*Par_mesh.grid_num(1),-1/3*Par_mesh.grid_num(1)], Par_mesh.grid_num(1),Par_mesh.grid_num(1));
        J3_noquasi = [ sparse(1/2*Par_mesh.grid_num(2)*Par_mesh.grid_num(1),1/2*Par_mesh.grid_num(2)*Par_mesh.grid_num(1)) , sparse(kron(speye(0.5*Par_mesh.grid_num(2)), J1_noquasi)) ;...
                       sparse(kron(speye(0.5*Par_mesh.grid_num(2)),J2_noquasi)), sparse(1/2*Par_mesh.grid_num(2)*Par_mesh.grid_num(1),1/2*Par_mesh.grid_num(2)*Par_mesh.grid_num(1)) ];
        
        Kx_forward_noquasi = spdiags(Ones_x,1,Par_mesh.grid_num(1),Par_mesh.grid_num(1));
        Kx_forward_noquasi(Par_mesh.grid_num(1),1) = exp(0);
        Kx_forward_noquasi = kron(Iz,kron(Iy,Kx_forward_noquasi));
        Ky_forward_noquasi = spdiags(Ones_y,1,Par_mesh.grid_num(2),Par_mesh.grid_num(2));
        Ky_forward_noquasi = kron(Ky_forward_noquasi,Ix);
        Ky_forward_noquasi((Par_mesh.grid_num(2)-1)*Par_mesh.grid_num(1)+1:end, 1:Par_mesh.grid_num(1)) = J1_noquasi;
        Ky_forward_noquasi = kron(Iz,Ky_forward_noquasi);
        Kz_forward_noquasi = spdiags(Ones_z,1,Par_mesh.grid_num(3),Par_mesh.grid_num(3));
        Kz_forward_noquasi = kron(Kz_forward_noquasi,kron(Iy,Ix));
        Kz_forward_noquasi((Par_mesh.grid_num(3)-1)*Par_mesh.grid_num(2)*Par_mesh.grid_num(1)+1:end, 1:Par_mesh.grid_num(2)*Par_mesh.grid_num(1)) = J3_noquasi;
                   
        B_eps_11_shift_x = Kx_forward_noquasi*B_eps_11;
        B_eps_12_shift_x = Kx_forward_noquasi*B_eps_12;
        B_eps_13_shift_x = Kx_forward_noquasi*B_eps_13;

        B_eps_21_shift_y = Ky_forward_noquasi*B_eps_21;
        B_eps_22_shift_y = Ky_forward_noquasi*B_eps_22;
        B_eps_23_shift_y = Ky_forward_noquasi*B_eps_23;

        B_eps_31_shift_z = Kz_forward_noquasi*B_eps_31;
        B_eps_32_shift_z = Kz_forward_noquasi*B_eps_32;
        B_eps_33_shift_z = Kz_forward_noquasi*B_eps_33;
end
B.B_eps_11 = spdiags( B_eps_11, 0, N, N);
B.B_eps_12 = spdiags( B_eps_12, 0, N, N);
B.B_eps_13 = spdiags( B_eps_13, 0, N, N);
B.B_eps_21 = spdiags( B_eps_21, 0, N, N);
B.B_eps_22 = spdiags( B_eps_22, 0, N, N);
B.B_eps_23 = spdiags( B_eps_23, 0, N, N);
B.B_eps_31 = spdiags( B_eps_31, 0, N, N);
B.B_eps_32 = spdiags( B_eps_32, 0, N, N);
B.B_eps_33 = spdiags( B_eps_33, 0, N, N);

B.B_eps_11_shift_x = spdiags( B_eps_11_shift_x, 0, N, N);
B.B_eps_12_shift_x = spdiags( B_eps_12_shift_x, 0, N, N);
B.B_eps_13_shift_x = spdiags( B_eps_13_shift_x, 0, N, N);

B.B_eps_21_shift_y = spdiags( B_eps_21_shift_y, 0, N, N);
B.B_eps_22_shift_y = spdiags( B_eps_22_shift_y, 0, N, N);
B.B_eps_23_shift_y = spdiags( B_eps_23_shift_y, 0, N, N);

B.B_eps_31_shift_z = spdiags( B_eps_31_shift_z, 0, N, N);
B.B_eps_32_shift_z = spdiags( B_eps_32_shift_z, 0, N, N);
B.B_eps_33_shift_z = spdiags( B_eps_33_shift_z, 0, N, N);
