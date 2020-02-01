function [ Popt ] = FAME_User_Option_gui( par )
%% Setting the parameters of domain mesh
% The grid numbers
Popt.mesh.x_grid_num = par.x_grid_num;
Popt.mesh.y_grid_num = par.x_grid_num;
Popt.mesh.z_grid_num = par.x_grid_num;

%% Setting the parameters of domain shape 
% The options for Popt.domain_shape are 'simple_cubic', 'face_centered_cubic', 'body_centered_cubic' or 'rhombohedral' 
Popt.lattice.lattice_type     = 'user_defined';
Popt.lattice.lattice_constant = par.lattice_constant; % do not change this term in 'user defined' case
Popt.material.sphere_radius   = par.sphere_radius;
Popt.material.cylinder_radius = par.cylinder_radius;
Popt.material.data_name       = par.file_name;
%% Setting the parameters of material
% The options for Popt.material.type are 'isotropic', 'anisotropic', 'biisotropic' or 'bianisotropic'
Popt.material.material_type = lower(par.mode);
% The permittivity and permeability of material
% When the material.type is 'isotropic' the constitutive relations should be given as 
%                                      D = £`*E,
%                                      B = £g*H.
% where £` is the eletric permittivity and £g is the magnetic permeability
switch Popt.material.material_type
    case {'isotropic','anisotropic'}
%     Popt.material.ele_permitt = [11 11 11 11];
%     Popt.material.mag_permeab = [1 1 1 1];
    Popt.material.ele_permitt = par.ele_permitt;
    Popt.material.mag_permeab = par.mag_permeab;
    
    case 'biisotropic'
    Popt.material.ele_permitt = 13;
    Popt.material.mag_permeab = 1;
    % When the material.type is 'biisotropic' or 'bianisotropic' the constitutive relations should be given as 
    %                                      D = £`*E + £i*H, 
    %                                      B = £g*H + £a*E.
    % where  £i = £q - i£e and £a = £q + i£e                                  
    %         and £q is the reciprocity of the medium and £e is the chirality of the medium
    Popt.material.reciprocity = par.reciprocity;
    Popt.material.chirality   = par.chirality;
    
end
%% Setting the parameters of graph
% This is the number of frequency (or eigenvalue) for each wave vector we want to plot on the graph,
Popt.graph.plot_num  = 8;

%% Setting the lattice in first Brilloion zone
% The number of paths in first Brillouin zone
% The partitions of each path which is connected by two wave vectors
Popt.lattice.part_num = par.lattice_part_num;
%% Setting the boundary condition
% The options for bd_cond of three directions are 'quasi_periodic', 'perfect_matching_layer', ( Prepariing... 'neumann' and 'dirichlet' ) 
Popt.boundary.bd_cond_x = 'quasi_periodic';
Popt.boundary.bd_cond_y = 'quasi_periodic';
Popt.boundary.bd_cond_z = 'quasi_periodic';

end