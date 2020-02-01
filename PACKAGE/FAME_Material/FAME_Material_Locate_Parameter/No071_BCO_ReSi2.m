function material = No071_BCO_ReSi2(sphere_radius,cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No071_BCO_ReSi2
%
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/AB2_oI6_71_a_i.html
%
% Edit at 2017/7/6 By Jia-Wei Lin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'ReSi2';
material.material_num = 2;

material.lattice_type = 'body_centered_orthorhombic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB2_oI6_71_a_i.html';

material.z2           = 0.339;

material.lattice_constant.a      = 3.144;
material.lattice_constant.b      = 0.994910941476*material.lattice_constant.a;
material.lattice_constant.c      = 2.44179389313*material.lattice_constant.a;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );

%% Shape description for Re (Sphere)
material.parameters{1}.name                 = 'Re';
material.parameters{1}.color_map            = [31 128 224]/255;
material.parameters{1}.sphere_centers       = [ 0, 0, 0 ];
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1));

%% Shape description for Si (Sphere)
material.parameters{2}.name                 = 'Si';
material.parameters{2}.color_map            = [221 169 104]/255;
material.parameters{2}.sphere_centers       = [   material.z2,   material.z2, 0;
                                                1-material.z2, 1-material.z2, 0];
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1));

%% Shape description for Ce  (Cylinders)
material.parameters{1}.cylinder_bot_centers = [ 0 0 0;
                                                0 0 0;
                                                0 0 0;
                                                1 1 1;
                                                1 1 1;
                                                1 1 1;
                                                0 1 1;
                                                0 1 1;
                                                1 0 0
                                                1 0 0];
material.parameters{1}.cylinder_top_centers = .5*[ [0 0 0] + material.parameters{2}.sphere_centers(1,:)
                                                   [0 0 0] + material.parameters{2}.sphere_centers(1,:) + [0 0 1]
                                                   [0 0 0] + material.parameters{2}.sphere_centers(2,:) + [0 0 1]
                                                   [1 1 1] + material.parameters{2}.sphere_centers(2,:) + [0 0 1]
                                                   [1 1 1] + material.parameters{2}.sphere_centers(1,:)
                                                   [1 1 1] + material.parameters{2}.sphere_centers(2,:)
                                                   [0 1 1] + material.parameters{2}.sphere_centers(1,:) + [0 0 1]
                                                   [0 1 1] + material.parameters{2}.sphere_centers(2,:) + [0 0 1]
                                                   [1 0 0] + material.parameters{2}.sphere_centers(1,:)
                                                   [1 0 0] + material.parameters{2}.sphere_centers(2,:)];
material.parameters{1}.cylinder_radius      = cylinder_radius(1)*ones( 1, size(material.parameters{1}.cylinder_bot_centers,1));
 
%% Shape description for Zn  (Cylinders)
material.parameters{2}.cylinder_bot_centers = [    material.parameters{2}.sphere_centers(1,:)
                                                   material.parameters{2}.sphere_centers(1,:) + [0 0 1]
                                                   material.parameters{2}.sphere_centers(2,:) + [0 0 1]
                                                   material.parameters{2}.sphere_centers(2,:) + [0 0 1]
                                                   material.parameters{2}.sphere_centers(1,:)
                                                   material.parameters{2}.sphere_centers(2,:)
                                                   material.parameters{2}.sphere_centers(1,:) + [0 0 1]
                                                   material.parameters{2}.sphere_centers(2,:) + [0 0 1]
                                                   material.parameters{2}.sphere_centers(1,:)
                                                   material.parameters{2}.sphere_centers(2,:)];
material.parameters{2}.cylinder_top_centers = .5*[ [0 0 0] + material.parameters{2}.sphere_centers(1,:)
                                                   [0 0 0] + material.parameters{2}.sphere_centers(1,:) + [0 0 1]
                                                   [0 0 0] + material.parameters{2}.sphere_centers(2,:) + [0 0 1]
                                                   [1 1 1] + material.parameters{2}.sphere_centers(2,:) + [0 0 1]
                                                   [1 1 1] + material.parameters{2}.sphere_centers(1,:)
                                                   [1 1 1] + material.parameters{2}.sphere_centers(2,:)
                                                   [0 1 1] + material.parameters{2}.sphere_centers(1,:) + [0 0 1]
                                                   [0 1 1] + material.parameters{2}.sphere_centers(2,:) + [0 0 1]
                                                   [1 0 0] + material.parameters{2}.sphere_centers(1,:)
                                                   [1 0 0] + material.parameters{2}.sphere_centers(2,:)];
material.parameters{2}.cylinder_radius      = cylinder_radius(2)*ones( 1, size(material.parameters{2}.cylinder_bot_centers,1));


%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};