function material = No180_Hexagonal_Mg2Ni(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No180_Hexagonal_Mg2Ni
%
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A2B_hP18_180_fi_bd.html
%
% Edit at 2017/7/25 by Jyun-Wei Lin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'Mg2Ni';
material.material_num = 2;

material.lattice_type = 'hexagonal';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A2B_hP18_180_fi_bd.html';

material.lattice_constant.a = 5.198;
material.lattice_constant.c = 2.54136206233*material.lattice_constant.a;
material.z3                 = 0.163;
material.x4                 = 0.1141;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );

%% Shape description for Ni
material.parameters{1}.name                 = 'Ni';
material.parameters{1}.color_map            = [0 255 0]/255;
material.parameters{1}.sphere_centers       = [          0,            0,               .5;
                                                         0,            0,              1/6;
                                                         0,            0,              5/6;
                                                        .5,            0,               .5;
                                                         0,           .5,              1/6;
                                                        .5,           .5,              5/6  ];
                                              
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 


%% Shape description for Mg
material.parameters{2}.name                 = 'Mg';
material.parameters{2}.color_map            = [163 255 70]/255;
material.parameters{2}.sphere_centers       = [            .5,               0,      material.z3;
                                                            0,              .5,  2/3+material.z3;
                                                           .5,              .5,  1/3+material.z3;
                                                           .5,               0,    1-material.z3;
                                                            0,              .5,  2/3-material.z3;
                                                           .5,              .5,  1/3-material.z3;
                                                  material.x4,   2*material.x4,                0;
                                              1-2*material.x4,   1-material.x4,              2/3;
                                                  material.x4,   1-material.x4,              1/3;
                                                1-material.x4, 1-2*material.x4,                0;
                                                2*material.x4,     material.x4,              2/3;
                                                1-material.x4,     material.x4,              1/3
                                              ];
%material.parameters{2}.sphere_centers=mod(material.parameters{2}.sphere_centers,1);
                                             
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 


%% Collect color map
material.color_map = {material.parameters{1}.color_map,material.parameters{2}.color_map};