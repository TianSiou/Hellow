function material = No063_CBCO_ZrSi2(sphere_radius,cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No063_CBCO_ZrSi2
%
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A2B_oC12_63_2c_c.html
%
% Edit at 2017/7/6 By Jia-Wei Lin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'ZrSi2';
material.material_num = 2;

material.lattice_type = 'c_base_centered_orthorhombic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A2B_oC12_63_2c_c.html';

material.y1           = 0.061;
material.y2           = 0.75;
material.y3           = 0.396;

material.lattice_constant.a      = 3.73;
material.lattice_constant.b      = 14.72;
material.lattice_constant.c      = 3.67;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );

%% Shape description for Ce (Sphere)
material.parameters{1}.name                 = 'Si';
material.parameters{1}.color_map            = [252 230 182]/255;
material.parameters{1}.sphere_centers       = [   1-material.y1,   material.y1, .25;
                                                    material.y1, 1-material.y1, .75;
                                                  1-material.y2,   material.y2, .25;
                                                    material.y2, 1-material.y2, .75];
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1));

%% Shape description for Zn (Sphere)
material.parameters{2}.name                 = 'Zn';
material.parameters{2}.color_map            = [41 233 104]/255;
material.parameters{2}.sphere_centers       = [   1-material.y3,   material.y3, .25;
                                                    material.y3, 1-material.y3, .75];
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1));

%% Shape description for Ce  (Cylinders)
material.parameters{1}.cylinder_bot_centers = [ material.parameters{1}.sphere_centers(1,:) - [1 0 0]
                                                material.parameters{1}.sphere_centers(1,:) + [0 1 0]
                                                material.parameters{1}.sphere_centers(4,:) - [1 0 0]
                                                material.parameters{1}.sphere_centers(4,:) + [0 1 0]
                                                material.parameters{1}.sphere_centers(3,:)
                                                material.parameters{1}.sphere_centers(2,:) - [0 1 0]
                                                material.parameters{1}.sphere_centers(2,:) + [1 0 0]
                                                material.parameters{1}.sphere_centers(4,:)
                                                material.parameters{1}.sphere_centers(1,:)
                                                material.parameters{1}.sphere_centers(2,:)
                                                material.parameters{1}.sphere_centers(2,:) - [0 1 0]
                                                material.parameters{1}.sphere_centers(2,:) + [1 0 0]
                                                material.parameters{1}.sphere_centers(3,:) - [0 1 0]
                                                material.parameters{1}.sphere_centers(3,:) + [1 0 0]
                                                material.parameters{1}.sphere_centers(1,:) + [0 1 0]
                                                material.parameters{1}.sphere_centers(1,:) - [1 0 0]];
material.parameters{1}.cylinder_top_centers = .5*[ material.parameters{1}.sphere_centers(1,:) + material.parameters{2}.sphere_centers(2,:) - [1 0 0]
                                                   material.parameters{1}.sphere_centers(1,:) + material.parameters{2}.sphere_centers(2,:) + [0 1 0]
                                                   material.parameters{1}.sphere_centers(4,:) + material.parameters{2}.sphere_centers(2,:) - [1 0 0]
                                                   material.parameters{1}.sphere_centers(4,:) + material.parameters{2}.sphere_centers(2,:) + [0 1 0]
                                                   material.parameters{1}.sphere_centers(3,:) + material.parameters{2}.sphere_centers(2,:)
                                                   material.parameters{1}.sphere_centers(2,:) + material.parameters{2}.sphere_centers(1,:) - [0 1 0]
                                                   material.parameters{1}.sphere_centers(2,:) + material.parameters{2}.sphere_centers(1,:) + [1 0 0]
                                                   material.parameters{1}.sphere_centers(4,:) + material.parameters{2}.sphere_centers(1,:)
                                                   material.parameters{1}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(1,:) - [-1 1 0]
                                                   material.parameters{1}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(1,:) + [-1 1 0]
                                                   material.parameters{1}.sphere_centers(2,:) + material.parameters{2}.sphere_centers(2,:) - [ 0 1 0]
                                                   material.parameters{1}.sphere_centers(2,:) + material.parameters{2}.sphere_centers(2,:) + [ 1 0 0]
                                                   material.parameters{1}.sphere_centers(3,:) + material.parameters{2}.sphere_centers(1,:) - [0 1 0]
                                                   material.parameters{1}.sphere_centers(3,:) + material.parameters{2}.sphere_centers(1,:) + [ 1 0 0]
                                                   material.parameters{1}.sphere_centers(1,:) + material.parameters{2}.sphere_centers(1,:) + [0 1 0]
                                                   material.parameters{1}.sphere_centers(1,:) + material.parameters{2}.sphere_centers(1,:) - [ 1 0 0]];
material.parameters{1}.cylinder_radius      = cylinder_radius(1)*ones( 1, size(material.parameters{1}.cylinder_bot_centers,1));
 
%% Shape description for Zn  (Cylinders)
material.parameters{2}.cylinder_bot_centers = [ material.parameters{2}.sphere_centers(2,:);
                                                material.parameters{2}.sphere_centers(2,:);
                                                material.parameters{2}.sphere_centers(2,:);
                                                material.parameters{2}.sphere_centers(2,:);
                                                material.parameters{2}.sphere_centers(2,:);
                                                material.parameters{2}.sphere_centers(1,:);
                                                material.parameters{2}.sphere_centers(1,:);
                                                material.parameters{2}.sphere_centers(1,:);
                                                material.parameters{2}.sphere_centers(2,:);
                                                material.parameters{2}.sphere_centers(2,:);
                                                material.parameters{2}.sphere_centers(1,:);
                                                material.parameters{2}.sphere_centers(1,:);
                                                material.parameters{2}.sphere_centers(1,:);
                                                material.parameters{2}.sphere_centers(1,:)];
material.parameters{2}.cylinder_top_centers = .5*[ material.parameters{1}.sphere_centers(1,:) + material.parameters{2}.sphere_centers(2,:) - [1 0 0]
                                                   material.parameters{1}.sphere_centers(1,:) + material.parameters{2}.sphere_centers(2,:) + [0 1 0]
                                                   material.parameters{1}.sphere_centers(4,:) + material.parameters{2}.sphere_centers(2,:) - [1 0 0]
                                                   material.parameters{1}.sphere_centers(4,:) + material.parameters{2}.sphere_centers(2,:) + [0 1 0]
                                                   material.parameters{1}.sphere_centers(3,:) + material.parameters{2}.sphere_centers(2,:)
                                                   material.parameters{1}.sphere_centers(2,:) + material.parameters{2}.sphere_centers(1,:) - [0 1 0]
                                                   material.parameters{1}.sphere_centers(2,:) + material.parameters{2}.sphere_centers(1,:) + [1 0 0]
                                                   material.parameters{1}.sphere_centers(4,:) + material.parameters{2}.sphere_centers(1,:)
                                                   material.parameters{1}.sphere_centers(2,:) + material.parameters{2}.sphere_centers(2,:) - [0 1 0]
                                                   material.parameters{1}.sphere_centers(2,:) + material.parameters{2}.sphere_centers(2,:) + [ 1 0 0]
                                                   material.parameters{1}.sphere_centers(3,:) + material.parameters{2}.sphere_centers(1,:) - [0 1 0]
                                                   material.parameters{1}.sphere_centers(3,:) + material.parameters{2}.sphere_centers(1,:) + [ 1 0 0]
                                                   material.parameters{1}.sphere_centers(1,:) + material.parameters{2}.sphere_centers(1,:) + [0 1 0]
                                                   material.parameters{1}.sphere_centers(1,:) + material.parameters{2}.sphere_centers(1,:) - [ 1 0 0]];
material.parameters{2}.cylinder_radius      = cylinder_radius(2)*ones( 1, size(material.parameters{2}.cylinder_bot_centers,1));


%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};