function material = No070_FCO_TiSi2(sphere_radius,cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No070_FCO_TiSi2
%
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A2B_oF24_70_e_a.html
%
% Edit at 2017/7/6 By Jia-Wei Lin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'TiSi2';
material.material_num = 2;

material.lattice_type = 'face_centered_orthorhombic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A2B_oF24_70_e_a.html';

material.x2           = 0.4615;

material.lattice_constant.a      = 8.2671;
material.lattice_constant.b      = 0.580614725841*material.lattice_constant.a;
material.lattice_constant.c      = 1.0342804611*material.lattice_constant.a;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Ti (Sphere)
material.parameters{1}.name                 = 'Ti';
material.parameters{1}.color_map            = [185 176 168]/255;
material.parameters{1}.sphere_centers       = [ .125, .125, .125 ;
                                                .875, .875, .875 ];
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1));

%% Shape description for Si (Sphere)
material.parameters{2}.name                 = 'Si';
material.parameters{2}.color_map            = [244 176 98]/255;
material.parameters{2}.sphere_centers       = [   1.25-material.x2,      material.x2,      material.x2;
                                                       material.x2, 1.25-material.x2, 1.25-material.x2;
                                                  -.25+material.x2,    1-material.x2,    1-material.x2;
                                                     1-material.x2, -.25+material.x2, -.25+material.x2];
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1));

%% Shape description for Ti  (Cylinders)
material.parameters{1}.cylinder_bot_centers = [ material.parameters{1}.sphere_centers(1,:) + [0 0 1]
                                                material.parameters{1}.sphere_centers(1,:) + [1 1 0]
                                                material.parameters{1}.sphere_centers(1,:) + [1 1 0]
                                                material.parameters{1}.sphere_centers(1,:)
                                                material.parameters{1}.sphere_centers(2,:)
                                                material.parameters{1}.sphere_centers(1,:)
                                                material.parameters{1}.sphere_centers(2,:)
                                                material.parameters{1}.sphere_centers(1,:)
                                                material.parameters{1}.sphere_centers(1,:) + [0 1 0]
                                                material.parameters{1}.sphere_centers(2,:) + [0 0 -1]
                                                material.parameters{1}.sphere_centers(1,:) + [0 0 1]
                                                material.parameters{1}.sphere_centers(2,:)
                                                material.parameters{1}.sphere_centers(1,:) + [0 0 1]
                                                material.parameters{1}.sphere_centers(2,:) + [0 -1 0]
                                                material.parameters{1}.sphere_centers(2,:) + [-1 -1 0]
                                                material.parameters{1}.sphere_centers(2,:) + [-1 -1 0]
                                                material.parameters{1}.sphere_centers(1,:) + [1 0 0]
                                                material.parameters{1}.sphere_centers(2,:) + [-1 0 0]
                                                material.parameters{1}.sphere_centers(2,:) + [0 0 -1]
                                                material.parameters{1}.sphere_centers(2,:) + [0 0 -1]];
material.parameters{1}.cylinder_top_centers = .5*[ material.parameters{2}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(1,:) + [0 0 1]
                                                   material.parameters{2}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(1,:) + [1 1 0]
                                                   material.parameters{2}.sphere_centers(1,:) + material.parameters{1}.sphere_centers(1,:) + [1 1 0]
                                                   material.parameters{2}.sphere_centers(4,:) + material.parameters{1}.sphere_centers(1,:)
                                                   material.parameters{2}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(2,:)
                                                   material.parameters{2}.sphere_centers(1,:) + material.parameters{1}.sphere_centers(1,:) + [-1 0  0]
                                                   material.parameters{2}.sphere_centers(3,:) + material.parameters{1}.sphere_centers(2,:) + [1 0 0]
                                                   material.parameters{2}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(1,:) + [0 -1 -1]
                                                   material.parameters{2}.sphere_centers(3,:) + material.parameters{1}.sphere_centers(1,:) + [0 1 0]
                                                   material.parameters{2}.sphere_centers(1,:) + material.parameters{1}.sphere_centers(2,:) + [0 0 -1]
                                                   material.parameters{2}.sphere_centers(3,:) + material.parameters{1}.sphere_centers(1,:) + [0 0 1]
                                                   material.parameters{1}.sphere_centers(2,:) + material.parameters{2}.sphere_centers(4,:) + [0 1 1]
                                                   material.parameters{2}.sphere_centers(1,:) + material.parameters{1}.sphere_centers(1,:) + [0 0 1]
                                                   material.parameters{2}.sphere_centers(1,:) + material.parameters{1}.sphere_centers(2,:) + [0 -1 0]
                                                   material.parameters{2}.sphere_centers(4,:) + material.parameters{1}.sphere_centers(2,:) + [-1 -1 0]
                                                   material.parameters{2}.sphere_centers(3,:) + material.parameters{1}.sphere_centers(2,:) + [-1 -1 0]
                                                   material.parameters{2}.sphere_centers(4,:) + material.parameters{1}.sphere_centers(1,:) + [1 0 0]
                                                   material.parameters{2}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(2,:) + [-1 0 0]
                                                   material.parameters{2}.sphere_centers(3,:) + material.parameters{1}.sphere_centers(2,:) + [0 0 -1]
                                                   material.parameters{2}.sphere_centers(4,:) + material.parameters{1}.sphere_centers(2,:) + [0 0 -1]];
material.parameters{1}.cylinder_radius      = cylinder_radius(1)*ones( 1, size(material.parameters{1}.cylinder_bot_centers,1));


%% Shape description for Si  (Cylinders)
material.parameters{2}.cylinder_bot_centers = [ material.parameters{2}.sphere_centers(1,:)
                                                material.parameters{2}.sphere_centers(1,:)
                                                material.parameters{2}.sphere_centers(1,:)
                                                material.parameters{2}.sphere_centers(2,:)
                                                material.parameters{2}.sphere_centers(2,:)
                                                material.parameters{2}.sphere_centers(2,:)
                                                material.parameters{2}.sphere_centers(2,:)
                                                material.parameters{2}.sphere_centers(1,:)
                                                material.parameters{2}.sphere_centers(4,:)
                                                material.parameters{2}.sphere_centers(2,:)
                                                material.parameters{2}.sphere_centers(3,:)
                                                material.parameters{2}.sphere_centers(2,:)
                                                material.parameters{2}.sphere_centers(4,:)
                                                material.parameters{2}.sphere_centers(2,:)
                                                material.parameters{2}.sphere_centers(4,:)
                                                material.parameters{2}.sphere_centers(3,:)
                                                material.parameters{2}.sphere_centers(4,:)
                                                material.parameters{2}.sphere_centers(3,:)
                                                material.parameters{2}.sphere_centers(4,:)
                                                material.parameters{2}.sphere_centers(4,:)
                                                material.parameters{2}.sphere_centers(3,:)
                                                material.parameters{2}.sphere_centers(4,:)
                                                material.parameters{2}.sphere_centers(2,:)
                                                material.parameters{2}.sphere_centers(1,:)
                                                material.parameters{2}.sphere_centers(3,:)
                                                material.parameters{2}.sphere_centers(1,:)
                                                material.parameters{2}.sphere_centers(1,:)
                                                material.parameters{2}.sphere_centers(3,:)
                                                material.parameters{2}.sphere_centers(1,:)
                                                material.parameters{2}.sphere_centers(2,:) + [0 -1 0]
                                                material.parameters{2}.sphere_centers(4,:)
                                                material.parameters{2}.sphere_centers(2,:)
                                                material.parameters{2}.sphere_centers(1,:)
                                                material.parameters{2}.sphere_centers(3,:)
                                                material.parameters{2}.sphere_centers(3,:)
                                                material.parameters{2}.sphere_centers(1,:) + [-1 0 0]
                                                material.parameters{2}.sphere_centers(3,:)
                                                material.parameters{2}.sphere_centers(4,:)
                                                material.parameters{2}.sphere_centers(3,:)
                                                material.parameters{2}.sphere_centers(4,:)];
material.parameters{2}.cylinder_top_centers = .5*[ material.parameters{2}.sphere_centers(1,:) + material.parameters{1}.sphere_centers(1,:) + [1 1 0]
                                                   material.parameters{2}.sphere_centers(1,:) + material.parameters{1}.sphere_centers(1,:) + [1 0  0]
                                                   material.parameters{2}.sphere_centers(1,:) + material.parameters{2}.sphere_centers(2,:)
                                                   material.parameters{2}.sphere_centers(1,:) + material.parameters{2}.sphere_centers(2,:)
                                                   material.parameters{2}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(1,:) + [0 0 1]
                                                   material.parameters{2}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(1,:) + [1 1 0]
                                                   material.parameters{2}.sphere_centers(2,:) + material.parameters{2}.sphere_centers(1,:) + [-1 0 1]
                                                   material.parameters{2}.sphere_centers(2,:) + material.parameters{2}.sphere_centers(1,:) + [1 0 -1]
                                                   
                                                   material.parameters{2}.sphere_centers(4,:) + material.parameters{1}.sphere_centers(1,:)
                                                   material.parameters{2}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(2,:)
                                                   
                                                   material.parameters{2}.sphere_centers(3,:) + material.parameters{1}.sphere_centers(2,:) + [ -1 0 0]
                                                   material.parameters{2}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(1,:) + [0 1 1]
                                                   material.parameters{2}.sphere_centers(2,:) + material.parameters{2}.sphere_centers(4,:) + [0 -1 0]
                                                   material.parameters{2}.sphere_centers(2,:) + material.parameters{2}.sphere_centers(4,:) + [0 1 0]
                                                   material.parameters{2}.sphere_centers(3,:) + material.parameters{2}.sphere_centers(4,:) + [1 0 -1]
                                                   material.parameters{2}.sphere_centers(3,:) + material.parameters{2}.sphere_centers(4,:) + [-1 0 1]
                                                   material.parameters{2}.sphere_centers(4,:) + material.parameters{2}.sphere_centers(3,:) 
                                                   material.parameters{2}.sphere_centers(3,:) + material.parameters{2}.sphere_centers(4,:)
                                                   material.parameters{2}.sphere_centers(4,:) + material.parameters{2}.sphere_centers(3,:) + [1 -1 -1]
                                                   material.parameters{2}.sphere_centers(4,:) + material.parameters{1}.sphere_centers(2,:) + [-1 -1 0]
                                                   material.parameters{2}.sphere_centers(3,:) + material.parameters{1}.sphere_centers(1,:) + [0 1 0]
                                                   material.parameters{2}.sphere_centers(4,:) + material.parameters{2}.sphere_centers(2,:) + [0 0 -1]
                                                   material.parameters{2}.sphere_centers(4,:) + material.parameters{2}.sphere_centers(2,:) + [0 0 1]
                                                   material.parameters{2}.sphere_centers(1,:) + material.parameters{1}.sphere_centers(2,:) + [0 0 -1]
                                                   material.parameters{2}.sphere_centers(3,:) + material.parameters{1}.sphere_centers(1,:) + [0 0 1]
                                                   material.parameters{2}.sphere_centers(1,:) + material.parameters{1}.sphere_centers(1,:) + [0 0 1]
                                                   material.parameters{2}.sphere_centers(1,:) + material.parameters{1}.sphere_centers(2,:) + [0 -1 0]
                                                   material.parameters{2}.sphere_centers(3,:) + material.parameters{1}.sphere_centers(2,:) + [-1 -1 0]
                                                   material.parameters{2}.sphere_centers(1,:) + material.parameters{2}.sphere_centers(2,:) + [0 -1 0]
                                                   material.parameters{2}.sphere_centers(1,:) + material.parameters{2}.sphere_centers(2,:) + [0 -1 0]
                                                   material.parameters{2}.sphere_centers(4,:) + material.parameters{1}.sphere_centers(1,:) + [1 0 0]
                                                   material.parameters{2}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(2,:) + [-1 0 0]
                                                   material.parameters{2}.sphere_centers(1,:) + material.parameters{2}.sphere_centers(3,:)
                                                   material.parameters{2}.sphere_centers(1,:) + material.parameters{2}.sphere_centers(3,:)
                                                   material.parameters{2}.sphere_centers(1,:) + material.parameters{2}.sphere_centers(3,:) + [-1 0 0]
                                                   material.parameters{2}.sphere_centers(1,:) + material.parameters{2}.sphere_centers(3,:) + [-1 0 0]
                                                   material.parameters{2}.sphere_centers(3,:) + material.parameters{1}.sphere_centers(2,:) + [0 0 -1]
                                                   material.parameters{2}.sphere_centers(4,:) + material.parameters{1}.sphere_centers(2,:) + [0 0 -1]
                                                   material.parameters{2}.sphere_centers(4,:) + material.parameters{2}.sphere_centers(3,:) + [0 1 0]
                                                   material.parameters{2}.sphere_centers(4,:) + material.parameters{2}.sphere_centers(3,:) + [0 -1 0]];
material.parameters{2}.cylinder_radius      = cylinder_radius(2)*ones( 1, size(material.parameters{2}.cylinder_bot_centers,1));


%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};