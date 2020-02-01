function [ Popt ] = FAME_User_Option( Popt )
%% Mesh setting
ncell = 27;
% n_supercell = [2 , 2];
Popt.mesh.grid_num = [8,8,8] .* [1,1,1] * (ncell.^(1/3) );
% Popt.mesh.grid_num = [16,16,16].*[sum(n_supercell), 1 , 1]; % The grid numbers
%  Popt.mesh.grid_num = [16,16,16].*[1,  1 , sum(n_supercell)];
Popt.eig.eigen_wanted  = 75;
%% Lattice settings
% -------------------- settings for 'user_defined' -----------------------
Popt.lattice.lattice_type     = 'user_defined';  % open for material
% ---- settings for 'FAME_Material_Isofunction_No230_BCC_DoubleGyroid' ---
% Popt.lattice.lattice_type       = 'body_centered_cubic';
% Popt.lattice.lattice_constant.a = 1;
% ----------- settings for 'FAME_Material_Isofunction_test' --------------
% Popt.lattice.lattice_type       = 'triclinic';
% SC_lattice_constant_a = 1;
% Popt.lattice.lattice_vec_a = [ 0.5*SC_lattice_constant_a*[-1;1;1], sum(n_supercell)*0.5*SC_lattice_constant_a*[1;-1;1], 0.5*SC_lattice_constant_a*[1;1;-1]];
% Popt.lattice.lattice_vec_a = [ sum(n_supercell)*0.5*SC_lattice_constant_a*[-1;1;1], sum(n_supercell)*0.5*SC_lattice_constant_a*[1;-1;1], sum(n_supercell)*0.5*SC_lattice_constant_a*[1;1;-1]];
%% Material settings
% ------------------ settings for 'user_defined' ------------------------
Popt.material.sphere_radius   = 0.3*[ 1 1 1 1 ];
Popt.material.cylinder_radius = 0.3*[ 1 1 1 1 ];
Popt.material.data_name       = 'No000_Primitive_Tetragonal_testDemo_cube'; % No227_FCC_DiamondStructure %%No221_SC_Demo%%cNo000_Primitive_Tetragonal_testDemo_cube
% ----------- settings for 'FAME_Material_Isofunction_test' --------------
% Popt.material.isofunction = @(x,y,z,a1,a2,a3) FAME_Material_Isofunction_BCC_SingleGyroid( x,y,z,a1,a2,a3,1 );
% Popt.material.isofunction = @(x,y,z,a1,a2,a3) FAME_Material_Isofunction_BCC_DoubleGyroid( x,y,z,a1,a2,a3,1 );
% Popt.material.isofunction = @(x,y,z,a1,a2,a3) FAME_Material_Isofunction_BCC_SingleGyroid_Defect( x,y,z,a1,a2,a3,1, 1);
% Popt.material.isofunction = @(x,y,z,a1,a2,a3) FAME_Material_Isofunction_BCC_DoubleGyroid_Defect( x,y,z,a1,a2,a3,1, 1);
% Popt.material.isofunction = @(x,y,z,a1,a2,a3) FAME_Material_Isofunction_TR_SGDG( x,y,z,a1,a2,a3,1, n_supercell);
% Popt.material.isofunction = @(x,y,z,a1,a2,a3) FAME_Material_Isofunction_TR_SGDG_Defect( x,y,z,a1,a2,a3,1, n_supercell );
Popt.material.isovalue    = [1.1, 1.1 -0.15];
Popt.material.periodic_judge   = 'on';
Popt.material.display_material = 'on';
Popt.material.material_type    = 'isotropic'; % The options for Popt.material.type are 'isotropic', 'anisotropic', 'biisotropic'

switch Popt.material.material_type
    case 'isotropic'
        % The permittivity and permeability of material
        % When the material.type is 'isotropic' the constitutive relations should be given as
        %                                      D = �`*E,
        %                                      B = �g*H.
        % where �` is the eletric permittivity and �g is the magnetic permeability
        Popt.material.ele_permitt_in    = [16 16 1];
        Popt.material.ele_permitt_out   = 1;
        Popt.material.mag_permeab_in    = [1 1 1];
        Popt.material.mag_permeab_out   = 1;
    case 'anisotropic'
        % The corresponding permittivity tensor under y-dircection external magnetic field
        extmag_intensity             = 0.875;
        permitt                      = 13;
        permitt_2                    = permitt*extmag_intensity;     %
        permitt_1                    = sqrt(permitt^2 + permitt_2^2);
        Popt.material.ele_permitt_in{1}  = [      permitt_1,               0, -1i*permitt_2;
                                                          0,         permitt,             0;
                                               1i*permitt_2,               0,     permitt_1];
        Popt.material.ele_permitt_out{1} = eye(3);
        Popt.material.mag_permeab_in{1}  = eye(3);
        Popt.material.mag_permeab_out{1} = eye(3);
    case 'biisotropic'
        Popt.material.ele_permitt_in  = 16; % epsilon
        Popt.material.mag_permeab_in  = 1;  % mu
        Popt.material.ele_permitt_out = 1; % epsilon
        Popt.material.mag_permeab_out = 1;
        % When the material.type is 'biisotropic' or 'bianisotropic' the constitutive relations should be given as
        %                                      D = �`*E + �i*H,
        %                                      B = �g*H + �a*E.
        % where  �i = �q - i�e and �a = �q + i�e
        %         and �q is the reciprocity of the medium and �e is the chirality of the medium
        Popt.material.reciprocity_in  = 0;
        Popt.material.chirality_in    = 1;
        Popt.material.reciprocity_out = 0;
        Popt.material.chirality_out   = 0;
end
%% Eigensolver settings
% Popt.eig.eigen_wanted  = 10;% This is the number of frequency (or eigenvalue) for each wave vector we want to plot on the graph,
%% Reciprocal lattice Settings
Popt.recip_lattice.part_num = 10; % The partitions of each path which is connected by two wave vectors

%% Boundary settings
% The options for bd_cond of three directions are 'quasi_periodic', ( Prepariing... 'perfect_matching_layer', 'neumann' and 'dirichlet' )
Popt.boundary.bd_cond_x = 'quasi_periodic';
Popt.boundary.bd_cond_y = 'quasi_periodic';
Popt.boundary.bd_cond_z = 'quasi_periodic';

end