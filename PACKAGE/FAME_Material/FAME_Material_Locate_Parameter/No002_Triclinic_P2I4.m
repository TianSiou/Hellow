function material = No002_Triclinic_P2I4(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No002_Triclinic_P2I4
%
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A2B_aP6_2_2i_i.html
%
% Edit at 2017/6/23 By Jia-Wei Lin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'P2I4';
material.material_num = 2;

material.lattice_type = 'triclinic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A2B_aP6_2_2i_i.html';
material.x1           = 0.557;
material.y1           = 0.73;
material.z1           = 0.165;

material.x2           = 0.82;
material.y2           = 0.803;
material.z2           = 0.695;

material.x3           = 0.397;
material.y3           = 0.639;
material.z3           = 0.463;

material.lattice_constant.a     = 4.56;
material.lattice_constant.b     = 7.06;
material.lattice_constant.c     = 7.4;
material.lattice_constant.alpha = pi*80.2/180;
material.lattice_constant.beta  = pi*106.96667/180;
material.lattice_constant.gamma = pi*98.2/180;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for C (Sphere)
material.parameters{1}.name                 = 'I';
material.parameters{1}.color_map            = [135 35 188]/255;
material.parameters{1}.sphere_centers       = [      material.x1,      material.y1,      material.z1;
                                                   1-material.x1,    1-material.y1,    1-material.z1;
                                                     material.x2,      material.y2,      material.z2;
                                                   1-material.x2,    1-material.y2,    1-material.z2 ];
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1));

%% Shape description for Fe (Sphere)
material.parameters{2}.name                 = 'P';
material.parameters{2}.color_map            = [213 90 49]/255;
material.parameters{2}.sphere_centers       = [      material.x3,      material.y3,      material.z3;
                                                   1-material.x3,    1-material.y3,    1-material.z3];
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1));

%% Shape description for Fe  (Cylinders)
material.parameters{2}.cylinder_bot_centers = [ material.parameters{2}.sphere_centers(1,:);
                                                material.parameters{2}.sphere_centers(2,:);
                                                material.parameters{2}.sphere_centers(1,:);
                                                material.parameters{2}.sphere_centers(2,:)];
material.parameters{2}.cylinder_top_centers = 0.5*( material.parameters{1}.sphere_centers + material.parameters{2}.cylinder_bot_centers );
material.parameters{2}.cylinder_radius      = cylinder_radius(2)*ones(1,size(material.parameters{2}.cylinder_bot_centers,1));

%% Shape description for C  (Cylinders)
material.parameters{1}.cylinder_bot_centers = material.parameters{1}.sphere_centers;
material.parameters{1}.cylinder_top_centers = 0.5*( material.parameters{1}.sphere_centers + material.parameters{2}.cylinder_bot_centers );
material.parameters{1}.cylinder_radius      = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_bot_centers,1));

%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};