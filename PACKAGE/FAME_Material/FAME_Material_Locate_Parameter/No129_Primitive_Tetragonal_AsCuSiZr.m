function material = No129_Primitive_Tetragonal_AsCuSiZr(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No129_Primitive_Tetragonal_AsCuSiZr
%
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/ABCD_tP8_129_c_b_a_c.html
%
% Edit at 2017/7/5 By Jia-Wei Lin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'AsCuSiZr';
material.material_num = 4;

material.lattice_type = 'primitive_tetragonal';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/ABCD_tP8_129_c_b_a_c.html';

material.z3           = 0.6793;
material.z4           = 0.2246;

material.lattice_constant.a      = 3.6736;
material.lattice_constant.c      = 2.60540069686*material.lattice_constant.a;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );

%% Shape description for Si (Sphere)
material.parameters{1}.name                 = 'Si';
material.parameters{1}.color_map            = [241 185 158]/255;
material.parameters{1}.sphere_centers       = [ .75, .25, 0;
                                                .25, .75, 0];
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1));

%% Shape description for Cu (Sphere)
material.parameters{2}.name                 = 'Cu';
material.parameters{2}.color_map            = [230 130 81]/255;
material.parameters{2}.sphere_centers       = [ .75, .25, .5;
                                                .25, .75, .5];
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1));

%% Shape description for As (Sphere)
material.parameters{3}.name                 = 'As';
material.parameters{3}.color_map            = [225 138 238]/255;
material.parameters{3}.sphere_centers       = [ .25, .25,   material.z3;
                                                .75, .75, 1-material.z3];
material.parameters{3}.sphere_radius        = sphere_radius(3)*ones(1,size(material.parameters{3}.sphere_centers,1));

%% Shape description for Zr (Sphere)
material.parameters{4}.name                 = 'Zr';
material.parameters{4}.color_map            = [115 225 222]/255;
material.parameters{4}.sphere_centers       = [ .25, .25,   material.z4;
                                                .75, .75, 1-material.z4];
                                            
material.parameters{4}.sphere_radius        = sphere_radius(4)*ones(1,size(material.parameters{4}.sphere_centers,1));

%% Shape description for Si  (Cylinders)
material.parameters{1}.cylinder_bot_centers = [ material.parameters{1}.sphere_centers(1,:);
                                                material.parameters{1}.sphere_centers(1,:);
                                                material.parameters{1}.sphere_centers(2,:);
                                                material.parameters{1}.sphere_centers(2,:);
                                                material.parameters{1}.sphere_centers(1,:)+[0 0 1];
                                                material.parameters{1}.sphere_centers(2,:)+[0 0 1]];
material.parameters{1}.cylinder_top_centers = .5*[ material.parameters{1}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(1,:);
                                                   material.parameters{4}.sphere_centers(1,:) + material.parameters{1}.sphere_centers(1,:);
                                                   material.parameters{4}.sphere_centers(1,:) + material.parameters{1}.sphere_centers(2,:);
                                                   material.parameters{1}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(1,:);
                                                   material.parameters{4}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(1,:) + [0 0 1];
                                                   material.parameters{4}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(2,:) + [0 0 1] ];
material.parameters{1}.cylinder_radius      = cylinder_radius(1)*ones( 1, size(material.parameters{1}.cylinder_bot_centers,1));

%% Shape description for Cu  (Cylinders)
material.parameters{2}.cylinder_bot_centers = [ material.parameters{2}.sphere_centers(1,:);
                                                material.parameters{2}.sphere_centers(1,:);
                                                material.parameters{2}.sphere_centers(2,:);
                                                material.parameters{2}.sphere_centers(2,:)];
material.parameters{2}.cylinder_top_centers = .5*[ material.parameters{2}.sphere_centers(1,:) + material.parameters{3}.sphere_centers(1,:);
                                                   material.parameters{2}.sphere_centers(1,:) + material.parameters{3}.sphere_centers(2,:);
                                                   material.parameters{2}.sphere_centers(2,:) + material.parameters{3}.sphere_centers(1,:);
                                                   material.parameters{2}.sphere_centers(2,:) + material.parameters{3}.sphere_centers(2,:)];
material.parameters{2}.cylinder_radius      = cylinder_radius(2)*ones( 1, size(material.parameters{2}.cylinder_bot_centers,1));


%% Shape description for As  (Cylinders)
material.parameters{3}.cylinder_bot_centers = [ material.parameters{3}.sphere_centers(1,:);
                                                material.parameters{3}.sphere_centers(2,:);
                                                material.parameters{3}.sphere_centers(1,:);
                                                material.parameters{3}.sphere_centers(2,:);
                                                material.parameters{3}.sphere_centers(1,:);
                                                material.parameters{3}.sphere_centers(2,:)];
material.parameters{3}.cylinder_top_centers = .5*[ material.parameters{4}.sphere_centers(2,:) + material.parameters{3}.sphere_centers(1,:)   ;
                                                   material.parameters{4}.sphere_centers(1,:) + material.parameters{3}.sphere_centers(2,:)   ;
                                                   material.parameters{2}.sphere_centers(1,:) + material.parameters{3}.sphere_centers(1,:);
                                                   material.parameters{2}.sphere_centers(1,:) + material.parameters{3}.sphere_centers(2,:);
                                                   material.parameters{2}.sphere_centers(2,:) + material.parameters{3}.sphere_centers(1,:);
                                                   material.parameters{2}.sphere_centers(2,:) + material.parameters{3}.sphere_centers(2,:)];
material.parameters{3}.cylinder_radius      = cylinder_radius(3)*ones( 1, size(material.parameters{3}.cylinder_bot_centers,1));


%% Shape description for Zr  (Cylinders)
material.parameters{4}.cylinder_bot_centers = [ material.parameters{4}.sphere_centers(1,:);
                                                material.parameters{4}.sphere_centers(1,:);
                                                material.parameters{4}.sphere_centers(2,:);
                                                material.parameters{4}.sphere_centers(2,:);
                                                material.parameters{4}.sphere_centers(1,:);
                                                material.parameters{4}.sphere_centers(2,:)];
material.parameters{4}.cylinder_top_centers = .5*[ material.parameters{4}.sphere_centers(1,:) + material.parameters{1}.sphere_centers(1,:)          ;
                                                   material.parameters{4}.sphere_centers(1,:) + material.parameters{1}.sphere_centers(2,:)          ;
                                                   material.parameters{4}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(1,:) + [0 0 1];
                                                   material.parameters{4}.sphere_centers(2,:) + material.parameters{1}.sphere_centers(2,:) + [0 0 1];
                                                   material.parameters{4}.sphere_centers(1,:) + material.parameters{3}.sphere_centers(2,:)   ;
                                                   material.parameters{4}.sphere_centers(2,:) + material.parameters{3}.sphere_centers(1,:)];
material.parameters{4}.cylinder_radius      = cylinder_radius(4)*ones( 1, size(material.parameters{4}.cylinder_bot_centers,1));


%% Collect color map
material.color_map = {material.parameters{1}.color_map,material.parameters{2}.color_map,material.parameters{3}.color_map,material.parameters{4}.color_map};