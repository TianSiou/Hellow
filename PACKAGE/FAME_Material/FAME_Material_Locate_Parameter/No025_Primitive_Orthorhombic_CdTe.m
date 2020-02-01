function material = No025_Primitive_Orthorhombic_CdTe(sphere_radius,cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No025_Primitive_Orthorhombic_CdTe
%
% This is a high-pressure phase of CdTe. We use the data given for a pressure of 19.3 GPa.
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/AB_oP2_25_b_a.html
%
% Edit at 2017/7/6 By Jia-Wei Lin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'CdTe';
material.material_num = 2;

material.lattice_type = 'primitive_orthorhombic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB_oP2_25_b_a.html';

material.z1           = 0.0;
material.z2           = 0.25;

material.lattice_constant.a      = 2.8102;
material.lattice_constant.b      = 1.87104120703*material.lattice_constant.a;
material.lattice_constant.c      = 1.0769696107*material.lattice_constant.a;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );

%% Shape description for Cd (Sphere)
material.parameters{1}.name                 = 'Cd';
material.parameters{1}.color_map            = [248 198 156]/255;
material.parameters{1}.sphere_centers       = [   0,   .5,   material.z2];
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1));

%% Shape description for Te (Sphere)
material.parameters{2}.name                 = 'Te';
material.parameters{2}.color_map            = [230 104 36]/255;
material.parameters{2}.sphere_centers       = [   0,   0,   material.z1];
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1));

%% Shape description for Cd  (Cylinders)
material.parameters{1}.cylinder_bot_centers = [   material.parameters{1}.sphere_centers;
                                                  material.parameters{1}.sphere_centers;
                                                  material.parameters{1}.sphere_centers;
                                                  material.parameters{1}.sphere_centers];
material.parameters{1}.cylinder_top_centers = .5*[ material.parameters{1}.sphere_centers + material.parameters{2}.sphere_centers
                                                   material.parameters{1}.sphere_centers + material.parameters{2}.sphere_centers + [0 1 0]
                                                   material.parameters{1}.sphere_centers + material.parameters{2}.sphere_centers + [0 0 1]
                                                   material.parameters{1}.sphere_centers + material.parameters{2}.sphere_centers + [0 1 1]];
material.parameters{1}.cylinder_radius      = cylinder_radius(1)*ones( 1, size(material.parameters{1}.cylinder_bot_centers,1));

%% Shape description for Te  (Cylinders)
material.parameters{2}.cylinder_bot_centers = [  material.parameters{2}.sphere_centers
                                                 material.parameters{2}.sphere_centers + [0 1 0]
                                                 material.parameters{2}.sphere_centers + [0 0 1]
                                                 material.parameters{2}.sphere_centers + [0 1 1]
                                                 material.parameters{2}.sphere_centers
                                                 material.parameters{2}.sphere_centers + [0 1 0]
                                                 material.parameters{2}.sphere_centers + [0 0 1]
                                                 material.parameters{2}.sphere_centers + [0 1 1]];
material.parameters{2}.cylinder_top_centers = .5*[ material.parameters{1}.sphere_centers + material.parameters{2}.sphere_centers
                                                   material.parameters{1}.sphere_centers + material.parameters{2}.sphere_centers + [0 1 0]
                                                   material.parameters{1}.sphere_centers + material.parameters{2}.sphere_centers + [0 0 1]
                                                   material.parameters{1}.sphere_centers + material.parameters{2}.sphere_centers + [0 1 1]
                                                   2*(material.parameters{2}.sphere_centers + [1 0 0])
                                                   2*(material.parameters{2}.sphere_centers + [1 1 0])
                                                   2*(material.parameters{2}.sphere_centers + [1 0 1])
                                                   2*(material.parameters{2}.sphere_centers + [1 1 1]) ];
material.parameters{2}.cylinder_radius      = cylinder_radius(2)*ones( 1, size(material.parameters{2}.cylinder_bot_centers,1));


%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};