function material = No139_BCT_Si(sphere_radius,cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No139_BCT_Si
%
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A_tI4_139_e.html
%
% Edit at 2017/7/5 By Jia-Wei Lin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'Si';
material.material_num = 1;

material.lattice_type = 'body_centered_tetragonal';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A_tI4_139_e.html';

material.z1           = 0.819;

material.lattice_constant.a      = 3.34916;
material.lattice_constant.c      = 1.94217355994*material.lattice_constant.a;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );

%% Shape description for Si (Sphere)
material.parameters{1}.name                 = 'Si';
material.parameters{1}.color_map            = [226 226 80]/255;
material.parameters{1}.sphere_centers       = [   0,   material.z1,   material.z1;
                                                  0, 1-material.z1, 1-material.z1 ];
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1));

%% Shape description for Si  (Cylinders)
material.parameters{1}.cylinder_bot_centers = [   material.parameters{1}.sphere_centers(2,:)
                                                  material.parameters{1}.sphere_centers(1,:) + [1 0 0]
                                                  material.parameters{1}.sphere_centers(1,:)
                                                  material.parameters{1}.sphere_centers(1,:) + [1 0 0]
                                                  material.parameters{1}.sphere_centers(1,:) + [1 0 0]
                                                  material.parameters{1}.sphere_centers(2,:)
                                                  material.parameters{1}.sphere_centers(1,:) + [1 0 0]
                                                  material.parameters{1}.sphere_centers(2,:)
                                                  material.parameters{1}.sphere_centers(2,:) + [1 0 0]
                                                  material.parameters{1}.sphere_centers(2,:) + [1 0 0]];
material.parameters{1}.cylinder_top_centers = .5*[ material.parameters{1}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(1,:) + [1 0 0];
                                                   material.parameters{1}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(1,:) + [1 0 0];
                                                   material.parameters{1}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(1,:) + [1 1 1];
                                                   material.parameters{1}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(1,:) + [2 1 0];
                                                   material.parameters{1}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(1,:) + [2 0 1];
                                                   material.parameters{1}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(1,:) - [0 0 1];
                                                   material.parameters{1}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(1,:) + [2 1 1];
                                                   material.parameters{1}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(1,:) - [0 1 1];
                                                   material.parameters{1}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(1,:) - [-1 1 1];
                                                   material.parameters{1}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(1,:) - [-2 1 0]];
material.parameters{1}.cylinder_radius      = cylinder_radius(1)*ones( 1, size(material.parameters{1}.cylinder_bot_centers,1));

%% Collect color map
material.color_map = {material.parameters{1}.color_map};