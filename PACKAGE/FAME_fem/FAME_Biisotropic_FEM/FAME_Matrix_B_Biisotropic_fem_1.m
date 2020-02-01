function B = FAME_Matrix_B_Biisotropic_fem_1( Par_mesh, Par_lattice, Par_material  )
%  (xi)£i = £q - i£e and (zeta)£a = £q + i£e

ele_permitt_in = Par_material.ele_permitt_in;
mag_permeab_in = Par_material.ele_permitt_out;
reciprocity_in = Par_material.reciprocity_in;
chirality_in   = Par_material.chirality_in;

N         = Par_mesh.grid_num(1)*Par_mesh.grid_num(2)*Par_mesh.grid_num(3);
inner_idx = Par_material.B.org_idx{1};

Ix = speye(Par_mesh.grid_num(1));
Iy = speye(Par_mesh.grid_num(2));
Iz = speye(Par_mesh.grid_num(3));
Ones_x = ones(Par_mesh.grid_num(1),1);
Ones_y = ones(Par_mesh.grid_num(2),1);
Ones_z = ones(Par_mesh.grid_num(3),1);

switch Par_lattice.lattice_type
    case 'simple_cubic'
        Kx_forward_p = spdiags(Ones_x,1,Par_mesh.grid_num(1),Par_mesh.grid_num(1));
        Kx_forward_p(Par_mesh.grid_num(1),1) = 1;
        Kx_forward_p = kron(Iz,kron(Iy,Kx_forward_p));
        Ky_forward_p = spdiags(Ones_y,1,Par_mesh.grid_num(2),Par_mesh.grid_num(2));
        Ky_forward_p(Par_mesh.grid_num(2),1) = 1;
        Ky_forward_p = kron(Iz,kron(Ky_forward_p,Ix));
        Kz_forward_p = spdiags(Ones_z,1,Par_mesh.grid_num(3),Par_mesh.grid_num(3));
        Kz_forward_p(Par_mesh.grid_num(3),1) = 1;
        Kz_forward_p = kron(Kz_forward_p,kron(Iy,Ix));

        Kxy_forward_p = Kx_forward_p*Ky_forward_p;
        Kyz_forward_p = Ky_forward_p*Kz_forward_p;
        Kzx_forward_p = Kz_forward_p*Kx_forward_p;

        B_eps = ones(N,1);    B_eps(inner_idx) = ele_permitt_in;
        B_mu   = ones(N,1);    B_mu(inner_idx)   = mag_permeab_in;
        B_zeta = zeros(N,1);   B_zeta(inner_idx) = reciprocity_in + 1i*chirality_in;
        B_xi   = zeros(N,1);   B_xi(inner_idx)   = reciprocity_in - 1i*chirality_in;

        B_eps_x  = Kx_forward_p*B_eps;   B_eps_y  = Ky_forward_p*B_eps;   B_eps_z  = Kz_forward_p*B_eps;
        % B_eps_xy = Kxy_forward_p*B_eps;  B_eps_yz = Kyz_forward_p*B_eps;  B_eps_zx = Kzx_forward_p*B_eps;
        B_mu_x    = Kx_forward_p*B_mu;     B_mu_y    = Ky_forward_p*B_mu;     B_mu_z    = Kz_forward_p*B_mu;
        B_mu_xy   = Kxy_forward_p*B_mu;    B_mu_yz   = Kyz_forward_p*B_mu;    B_mu_zx   = Kzx_forward_p*B_mu;
        B_zeta_x  = Kx_forward_p*B_zeta;   B_zeta_y  = Ky_forward_p*B_zeta;   B_zeta_z  = Kz_forward_p*B_zeta;
        B_zeta_xy = Kxy_forward_p*B_zeta;  B_zeta_yz = Kyz_forward_p*B_zeta;  B_zeta_zx = Kzx_forward_p*B_zeta;
        B_xi_x    = Kx_forward_p*B_xi;     B_xi_y    = Ky_forward_p*B_xi;     B_xi_z    = Kz_forward_p*B_xi;
        % B_xi_xy   = Kxy_forward_p*B_xi;    B_xi_yz   = Kyz_forward_p*B_xi;    B_xi_zx   = Kzx_forward_p*B_xi;

        B.B_eps = spdiags(B_eps,0,N,N);
        B.B_mu = spdiags(B_mu,0,N,N);
        B.B_zeta = spdiags(B_zeta,0,N,N);
        B.B_xi = spdiags(B_xi,0,N,N);
        B.B_eps_x  = spdiags(B_eps_x,0,N,N);   B.B_eps_y  = spdiags(B_eps_y,0,N,N);    B.B_eps_z  = spdiags(B_eps_z,0,N,N);
        % B_eps_xy = spdiags(B_eps_xy,0,N,N);  B_eps_yz = spdiags(B_eps_yz,0,N,N);   B_eps_zx = spdiags(B_eps_zx,0,N,N);
        B.B_mu_x    = spdiags(B_mu_x,0,N,N);     B.B_mu_y    = spdiags(B_mu_y,0,N,N);      B.B_mu_z    = spdiags(B_mu_z,0,N,N);
        B.B_mu_xy   = spdiags(B_mu_xy,0,N,N);    B.B_mu_yz   = spdiags(B_mu_yz,0,N,N);     B.B_mu_zx   = spdiags(B_mu_zx,0,N,N);
        B.B_zeta_x  = spdiags(B_zeta_x,0,N,N);   B.B_zeta_y  = spdiags(B_zeta_y,0,N,N);    B.B_zeta_z  = spdiags(B_zeta_z,0,N,N);
        B.B_zeta_xy = spdiags(B_zeta_xy,0,N,N);  B.B_zeta_yz = spdiags(B_zeta_yz,0,N,N);   B.B_zeta_zx = spdiags(B_zeta_zx,0,N,N);
        B.B_xi_x    = spdiags(B_xi_x,0,N,N);     B.B_xi_y    = spdiags(B_xi_y,0,N,N);      B.B_xi_z    = spdiags(B_xi_z,0,N,N);
        % B_xi_xy   = spdiags(B_xi_xy,0,N,N);    B_xi_yz   = spdiags(B_xi_yz,0,N,N);     B_xi_zx   = spdiags(B_xi_zx,0,N,N);
    case 'face_centered_cubic'
        J2_p =   spdiags([ Ones_x, exp(0)*Ones_x ], [-0.5*Par_mesh.grid_num(1),0.5*Par_mesh.grid_num(1)], Par_mesh.grid_num(1),Par_mesh.grid_num(1));
        J3_p = [ sparse((Par_mesh.grid_num(2)*Par_mesh.grid_num(1))/3,2*(Par_mesh.grid_num(2)*Par_mesh.grid_num(1)/3)),    exp(0)*speye((Par_mesh.grid_num(2)*Par_mesh.grid_num(1)/3)) ;...
                       kron(speye(2*Par_mesh.grid_num(2)/3),J2_p)                                     ,    sparse(2*(Par_mesh.grid_num(2)*Par_mesh.grid_num(1)/3),(Par_mesh.grid_num(2)*Par_mesh.grid_num(1))/3) ];
        Kx_forward_p = spdiags(Ones_x,1,Par_mesh.grid_num(1),Par_mesh.grid_num(1));
        Kx_forward_p(Par_mesh.grid_num(1),1) = exp(0);
        Kx_forward_p = kron(Iz,kron(Iy,Kx_forward_p));
        Ky_forward_p = spdiags(Ones_y,1,Par_mesh.grid_num(2),Par_mesh.grid_num(2));
        Ky_forward_p = kron(Ky_forward_p,Ix);
        Ky_forward_p((Par_mesh.grid_num(2)-1)*Par_mesh.grid_num(1)+1:end, 1:Par_mesh.grid_num(1)) = J2_p;
        Ky_forward_p = kron(Iz,Ky_forward_p);
        Kz_forward_p = spdiags(Ones_z,1,Par_mesh.grid_num(3),Par_mesh.grid_num(3));
        Kz_forward_p = kron(Kz_forward_p,kron(Iy,Ix));
        Kz_forward_p((Par_mesh.grid_num(3)-1)*Par_mesh.grid_num(2)*Par_mesh.grid_num(1)+1:end, 1:Par_mesh.grid_num(2)*Par_mesh.grid_num(1)) = J3_p;
         
        Kxy_forward_p = Kx_forward_p*Ky_forward_p;
        Kyz_forward_p = Ky_forward_p*Kz_forward_p;
        Kzx_forward_p = Kz_forward_p*Kx_forward_p;

        B_eps = ones(N,1);    B_eps(inner_idx) = ele_permitt_in;
        B_mu   = ones(N,1);    B_mu(inner_idx)   = mag_permeab_in;
        B_zeta = zeros(N,1);   B_zeta(inner_idx) = reciprocity_in + 1i*chirality_in;
        B_xi   = zeros(N,1);   B_xi(inner_idx)   = reciprocity_in - 1i*chirality_in;

        B_eps_x  = Kx_forward_p*B_eps;   B_eps_y  = Ky_forward_p*B_eps;   B_eps_z  = Kz_forward_p*B_eps;
        % B_eps_xy = Kxy_forward_p*B_eps;  B_eps_yz = Kyz_forward_p*B_eps;  B_eps_zx = Kzx_forward_p*B_eps;
        B_mu_x    = Kx_forward_p*B_mu;     B_mu_y    = Ky_forward_p*B_mu;     B_mu_z    = Kz_forward_p*B_mu;
        B_mu_xy   = Kxy_forward_p*B_mu;    B_mu_yz   = Kyz_forward_p*B_mu;    B_mu_zx   = Kzx_forward_p*B_mu;
        B_zeta_x  = Kx_forward_p*B_zeta;   B_zeta_y  = Ky_forward_p*B_zeta;   B_zeta_z  = Kz_forward_p*B_zeta;
        B_zeta_xy = Kxy_forward_p*B_zeta;  B_zeta_yz = Kyz_forward_p*B_zeta;  B_zeta_zx = Kzx_forward_p*B_zeta;
        B_xi_x    = Kx_forward_p*B_xi;     B_xi_y    = Ky_forward_p*B_xi;     B_xi_z    = Kz_forward_p*B_xi;
        % B_xi_xy   = Kxy_forward_p*B_xi;    B_xi_yz   = Kyz_forward_p*B_xi;    B_xi_zx   = Kzx_forward_p*B_xi;

        B.B_eps = spdiags(B_eps,0,N,N);
        B.B_mu = spdiags(B_mu,0,N,N);
        B.B_zeta = spdiags(B_zeta,0,N,N);
        B.B_xi = spdiags(B_xi,0,N,N);
        B.B_eps_x  = spdiags(B_eps_x,0,N,N);   B.B_eps_y  = spdiags(B_eps_y,0,N,N);    B.B_eps_z  = spdiags(B_eps_z,0,N,N);
        % B_eps_xy = spdiags(B_eps_xy,0,N,N);  B_eps_yz = spdiags(B_eps_yz,0,N,N);   B_eps_zx = spdiags(B_eps_zx,0,N,N);
        B.B_mu_x    = spdiags(B_mu_x,0,N,N);     B.B_mu_y    = spdiags(B_mu_y,0,N,N);      B.B_mu_z    = spdiags(B_mu_z,0,N,N);
        B.B_mu_xy   = spdiags(B_mu_xy,0,N,N);    B.B_mu_yz   = spdiags(B_mu_yz,0,N,N);     B.B_mu_zx   = spdiags(B_mu_zx,0,N,N);
        B.B_zeta_x  = spdiags(B_zeta_x,0,N,N);   B.B_zeta_y  = spdiags(B_zeta_y,0,N,N);    B.B_zeta_z  = spdiags(B_zeta_z,0,N,N);
        B.B_zeta_xy = spdiags(B_zeta_xy,0,N,N);  B.B_zeta_yz = spdiags(B_zeta_yz,0,N,N);   B.B_zeta_zx = spdiags(B_zeta_zx,0,N,N);
        B.B_xi_x    = spdiags(B_xi_x,0,N,N);     B.B_xi_y    = spdiags(B_xi_y,0,N,N);      B.B_xi_z    = spdiags(B_xi_z,0,N,N);
        % B_xi_xy   = spdiags(B_xi_xy,0,N,N);    B_xi_yz   = spdiags(B_xi_yz,0,N,N);     B_xi_zx   = spdiags(B_xi_zx,0,N,N);
end

end