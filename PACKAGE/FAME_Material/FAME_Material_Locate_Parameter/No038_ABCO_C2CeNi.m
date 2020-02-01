function material = No038_ABCO_C2CeNi(sphere_radius,cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No038_ABCO_C2CeNi
%
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A2BC_oC8_38_e_a_b.html
%
% Edit at 2017/7/6 By Jia-Wei Lin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'C2CeNi';
material.material_num = 3;

material.lattice_type = 'a_base_centered_orthorhombic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A2BC_oC8_38_e_a_b.html';

material.z1           = 0.0;
material.z2           = 0.6144;
material.y3           = 0.155;
material.z3           = 0.2914;

material.lattice_constant.a      = 3.875;
material.lattice_constant.b      = 1.17470967742*material.lattice_constant.a;
material.lattice_constant.c      = 1.59019354839*material.lattice_constant.a;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Ce (Sphere)
material.parameters{1}.name                 = 'Ce';
material.parameters{1}.color_map            = [252 230 182]/255;
material.parameters{1}.sphere_centers       = [   0, -material.z1,   material.z1];
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1));

%% Shape description for Ni (Sphere)
material.parameters{2}.name                 = 'Ni';
material.parameters{2}.color_map            = [41 233 104]/255;
material.parameters{2}.sphere_centers       = [   .5,   1-material.z2,   material.z2];
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1));

%% Shape description for C (Sphere)
material.parameters{3}.name                 = 'C';
material.parameters{3}.color_map            = [196 189 176]/255;
material.parameters{3}.sphere_centers       = [   .5,   1+material.y3-material.z3,   material.y3+material.z3
                                                  .5,   1-material.y3-material.z3,   material.z3-material.y3];
material.parameters{3}.sphere_radius        = sphere_radius(3)*ones(1,size(material.parameters{3}.sphere_centers,1));

%% Shape description for Ce  (Cylinders)
material.parameters{1}.cylinder_bot_centers = [   material.parameters{1}.sphere_centers + [0 1 0];
                                                  material.parameters{1}.sphere_centers + [0 1 1];
                                                  material.parameters{1}.sphere_centers + [1 1 0];
                                                  material.parameters{1}.sphere_centers + [1 1 1];
                                                  material.parameters{1}.sphere_centers;
                                                  material.parameters{1}.sphere_centers + [1 0 0];
                                                  material.parameters{1}.sphere_centers + [0 1 0];
                                                  material.parameters{1}.sphere_centers + [1 1 0]];
material.parameters{1}.cylinder_top_centers = .5*[ material.parameters{3}.sphere_centers(1,:) + material.parameters{1}.sphere_centers + [0 1 0];
                                                   material.parameters{3}.sphere_centers(1,:) + material.parameters{1}.sphere_centers + [0 1 1];
                                                   material.parameters{3}.sphere_centers(1,:) + material.parameters{1}.sphere_centers + [1 1 0];
                                                   material.parameters{3}.sphere_centers(1,:) + material.parameters{1}.sphere_centers + [1 1 1];
                                                   material.parameters{3}.sphere_centers(2,:) + material.parameters{1}.sphere_centers          ;
                                                   material.parameters{3}.sphere_centers(2,:) + material.parameters{1}.sphere_centers + [1 0 0];
                                                   material.parameters{3}.sphere_centers(2,:) + material.parameters{1}.sphere_centers + [0 1 0];
                                                   material.parameters{3}.sphere_centers(2,:) + material.parameters{1}.sphere_centers + [1 1 0]];
material.parameters{1}.cylinder_radius      = cylinder_radius(1)*ones( 1, size(material.parameters{1}.cylinder_bot_centers,1));


%% Shape description for Ni  (Cylinders)
material.parameters{2}.cylinder_bot_centers = [   material.parameters{2}.sphere_centers;
                                                  material.parameters{2}.sphere_centers;
                                                  material.parameters{2}.sphere_centers;
                                                  material.parameters{2}.sphere_centers];
material.parameters{2}.cylinder_top_centers = .5*[ material.parameters{3}.sphere_centers(1,:) + material.parameters{2}.sphere_centers - [0 1 0];
                                                   material.parameters{3}.sphere_centers(2,:) + material.parameters{2}.sphere_centers + [0 0 1];
                                                   material.parameters{3}.sphere_centers(1,:) + material.parameters{2}.sphere_centers;
                                                   material.parameters{3}.sphere_centers(2,:) + material.parameters{2}.sphere_centers];
material.parameters{2}.cylinder_radius      = cylinder_radius(2)*ones( 1, size(material.parameters{2}.cylinder_bot_centers,1));

%% Shape description for C  (Cylinders)
material.parameters{3}.cylinder_bot_centers = [   material.parameters{3}.sphere_centers(1,:);
                                                  material.parameters{3}.sphere_centers(2,:);
                                                  material.parameters{3}.sphere_centers(1,:);
                                                  material.parameters{3}.sphere_centers(2,:);
                                                  material.parameters{3}.sphere_centers(1,:) - [0 1 0];
                                                  material.parameters{3}.sphere_centers(2,:) + [0 0 1];
                                                  material.parameters{3}.sphere_centers(1,:);
                                                  material.parameters{3}.sphere_centers(1,:);
                                                  material.parameters{3}.sphere_centers(1,:);
                                                  material.parameters{3}.sphere_centers(1,:);
                                                  material.parameters{3}.sphere_centers(2,:);
                                                  material.parameters{3}.sphere_centers(2,:);
                                                  material.parameters{3}.sphere_centers(2,:);
                                                  material.parameters{3}.sphere_centers(2,:)];
material.parameters{3}.cylinder_top_centers = .5*[ material.parameters{3}.sphere_centers(1,:) + material.parameters{2}.sphere_centers;
                                                   material.parameters{3}.sphere_centers(2,:) + material.parameters{2}.sphere_centers;
                                                   material.parameters{3}.sphere_centers(1,:) + material.parameters{3}.sphere_centers(2,:);
                                                   material.parameters{3}.sphere_centers(2,:) + material.parameters{3}.sphere_centers(1,:);
                                                   material.parameters{3}.sphere_centers(1,:) + material.parameters{2}.sphere_centers - [0 1 0];
                                                   material.parameters{3}.sphere_centers(2,:) + material.parameters{2}.sphere_centers + [0 0 1];
                                                   material.parameters{3}.sphere_centers(1,:) + material.parameters{1}.sphere_centers + [0 1 0];
                                                   material.parameters{3}.sphere_centers(1,:) + material.parameters{1}.sphere_centers + [0 1 1];
                                                   material.parameters{3}.sphere_centers(1,:) + material.parameters{1}.sphere_centers + [1 1 0];
                                                   material.parameters{3}.sphere_centers(1,:) + material.parameters{1}.sphere_centers + [1 1 1];
                                                   material.parameters{3}.sphere_centers(2,:) + material.parameters{1}.sphere_centers          ;
                                                   material.parameters{3}.sphere_centers(2,:) + material.parameters{1}.sphere_centers + [1 0 0];
                                                   material.parameters{3}.sphere_centers(2,:) + material.parameters{1}.sphere_centers + [0 1 0];
                                                   material.parameters{3}.sphere_centers(2,:) + material.parameters{1}.sphere_centers + [1 1 0]];
material.parameters{3}.cylinder_radius      = cylinder_radius(3)*ones( 1, size(material.parameters{3}.cylinder_bot_centers,1));


%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map, material.parameters{3}.color_map};