function material = No012_BCM_AuTe2(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No012_BCM_AuTe2
%
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/AB2_mC6_12_a_i.html
%
% Edit at 2017/7/5 By Jia-Wei Lin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'AuTe2';
material.material_num = 2;

material.lattice_type = 'base_centered_monoclinic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB2_mC6_12_a_i.html';

material.x2           = 0.6879;
material.z2           = 0.2889;

material.lattice_constant.a      = 7.189;
material.lattice_constant.b      = 0.705105021561*material.lattice_constant.a;
material.lattice_constant.c      = 0.613019891501*material.lattice_constant.a;
material.lattice_constant.gamma  = pi*90.04/180;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );

%% Shape description for Au (Sphere)
material.parameters{1}.name                 = 'Au';
material.parameters{1}.color_map            = [226 226 80]/255;
material.parameters{1}.sphere_centers       = [ 0, 0, 0 ];
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1));

%% Shape description for Te (Sphere)
material.parameters{2}.name                 = 'Te';
material.parameters{2}.color_map            = [230 130 81]/255;
material.parameters{2}.sphere_centers       = [   material.x2,   material.z2,   material.x2;
                                                1-material.x2, 1-material.z2, 1-material.x2];
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1));

%% Shape description for Au  (Cylinders)
material.parameters{1}.cylinder_bot_centers = [ 0,  1,  0;
                                                1,  0,  1;
                                                1,  1,  1;
                                                0,  1,  1;
                                                0,  0,  0;
                                                1,  0,  0];
material.parameters{1}.cylinder_top_centers = .5*[ material.parameters{2}.sphere_centers(2,:) + [0 1 0];
                                                   material.parameters{2}.sphere_centers(1,:) + [1 0 1];
                                                   material.parameters{2}.sphere_centers(1,:) + [1 1 1];
                                                   material.parameters{2}.sphere_centers(2,:) + [0 1 1];
                                                   material.parameters{2}.sphere_centers(2,:) + [0 0 0];
                                                   material.parameters{2}.sphere_centers(1,:) + [1 0 0]];
material.parameters{1}.cylinder_radius      = cylinder_radius(1)*ones( 1, size(material.parameters{1}.cylinder_bot_centers,1));

%% Shape description for Te  (Cylinders)
material.parameters{2}.cylinder_bot_centers = [ material.parameters{2}.sphere_centers(2,:);
                                                material.parameters{2}.sphere_centers(1,:);
                                                material.parameters{2}.sphere_centers(1,:);
                                                material.parameters{2}.sphere_centers(2,:);
                                                material.parameters{2}.sphere_centers(2,:);
                                                material.parameters{2}.sphere_centers(1,:)];
material.parameters{2}.cylinder_top_centers = [ .5*(material.parameters{2}.sphere_centers(2,:) + [0 1 0]);
                                                .5*(material.parameters{2}.sphere_centers(1,:) + [1 0 1]);
                                                .5*(material.parameters{2}.sphere_centers(1,:) + [1 1 1]);
                                                .5*(material.parameters{2}.sphere_centers(2,:) + [0 1 1]);
                                                .5*(material.parameters{2}.sphere_centers(2,:) + [0 0 0]);
                                                .5*(material.parameters{2}.sphere_centers(1,:) + [1 0 0])];
material.parameters{2}.cylinder_radius      = cylinder_radius(2)*ones( 1, size(material.parameters{2}.cylinder_bot_centers,1));


%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};