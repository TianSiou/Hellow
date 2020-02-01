function material = No004_Primitive_Monclinic_Te(sphere_radius,cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No004_Primitive_Monclinic_Te
%
% This is a high-pressure phase of Te, stable in the 4-7 GPa range. The ground state of 
% Te appears to be £^¡VSe (A8).
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A_mP4_4_2a.html
%
% Edit at 2017/7/5 By Jia-Wei Lin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'Te';
material.material_num = 1;
material.lattice_type = 'primitive_monoclinic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A_mP4_4_2a.html';

material.x1           = 0.25;
material.y1           = 0.23;
material.z1           = 0.48;
material.x2           = 0.48;
material.y2           = 0.0;
material.z2           = 0.02;

material.lattice_constant.a      = 4.76;
material.lattice_constant.b      = 3.104;
material.lattice_constant.c      = 7.513;
material.lattice_constant.alpha  = pi*(180-92.71)/180;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );

permutation = [1 2 3];

%% Shape description for Te (Sphere)
material.parameters{1}.name                 = 'Te';
material.parameters{1}.color_map            = [135 35 188]/255;
material.parameters{1}.sphere_centers       = [      material.z1,      material.x1,      material.y1;
                                                   1-material.z1,   .5+material.x1,    1-material.y1;
                                                     material.z2,      material.x2,      material.y2;
                                                   1-material.z2,    1-material.x2,   .5+material.y2];

material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1));
material.parameters{1}.sphere_centers       = material.parameters{1}.sphere_centers(:,permutation);

%% Shape description for Te  (Cylinders)
material.parameters{1}.cylinder_bot_centers = [ material.parameters{1}.sphere_centers(3,:);
                                                material.parameters{1}.sphere_centers(1,:);
                                                material.parameters{1}.sphere_centers(4,:);
                                                material.parameters{1}.sphere_centers(2,:)];
material.parameters{1}.cylinder_top_centers = [ material.parameters{1}.sphere_centers(1,:);
                                                material.parameters{1}.sphere_centers(4,:);
                                                material.parameters{1}.sphere_centers(2,:);
                                                material.parameters{1}.sphere_centers(3,:) + [0 0 1] ];
material.parameters{1}.cylinder_radius      = cylinder_radius(1)*ones( 1, size(material.parameters{1}.cylinder_bot_centers,1));
material.parameters{1}.cylinder_bot_centers = material.parameters{1}.cylinder_bot_centers(:,permutation);
material.parameters{1}.cylinder_top_centers = material.parameters{1}.cylinder_top_centers(:,permutation);
%% Collect color map
material.color_map = {material.parameters{1}.color_map};