function material = No014_Primitive_Monoclinic_ZrO2(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No014_Primitive_Monoclinic_ZrO2
%
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A2B_mP12_14_2e_e.html
%
% Edit at 2017/8/16 By Hsiao-Han Huang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
material.Prototype    = 'ZrO2';
material.material_num = 2;

material.lattice_type = 'primitive_monoclinic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A2B_mP12_14_2e_e.html';
material.x1           = 0.07;
material.y1           = 0.3317;
material.z1           = 0.3447;

material.x2           = 0.4496;
material.y2           = 0.7569;
material.z2           = 0.4792;

material.x3           = 0.2754;
material.y3           = 0.0395;
material.z3           = 0.2083;

material.lattice_constant.b = 5.1505;
material.lattice_constant.a = 1.01186292593*material.lattice_constant.b;
material.lattice_constant.c = 1.03238520532*material.lattice_constant.b;
material.lattice_constant.alpha = pi*(180-99.23)/180;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );

%% Shape description for ¢Ý (Sphere)
material.parameters{1}.name                 = 'O';
material.parameters{1}.color_map            = [255 0 0]/255;

material.parameters{1}.sphere_centers       = [        material.x1,          material.y1,        material.z1;
                                                    -material.x1+1,       .5+material.y1,     .5-material.z1;
                                                    -material.x1+1,       -material.y1+1,     -material.z1+1;
                                                       material.x1,       .5-material.y1,     .5+material.z1;
                                                       material.x2,          material.y2,        material.z2;
                                                    -material.x2+1,     .5+material.y2-1,     .5-material.z2;
                                                    -material.x2+1,       -material.y2+1,     -material.z2+1;
                                                       material.x2,     .5-material.y2+1,     .5+material.z2 ];
material.parameters{1}.sphere_centers = material.parameters{1}.sphere_centers(:,[2,1,3]);                                                   
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(size(material.parameters{1}.sphere_centers,1),1);                                              

%% Shape description for Zr (Sphere)
material.parameters{2}.name                 = 'Zr';
material.parameters{2}.color_map            = [128 255 255]/255;

material.parameters{2}.sphere_centers       = [       material.x3,        material.y3,        material.z3;
                                                   -material.x3+1,     .5+material.y3,     .5-material.z3;
                                                   -material.x3+1,     -material.y3+1,     -material.z3+1;
                                                      material.x3,     .5-material.y3,     .5+material.z3;
                                              ];
material.parameters{2}.sphere_centers = material.parameters{2}.sphere_centers(:,[2,1,3]);                                          
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(size(material.parameters{2}.sphere_centers,1),1);                                              

%% Shape description for O  (Cylinders)
material.parameters{1}.cylinder_bot_centers = [        material.x1,          material.y1,        material.z1;
                                                    -material.x2+1,     .5+material.y2-1,     .5-material.z2;
                                                    -material.x2+1,       -material.y2+1,     -material.z2+1;
                                                    -material.x1+1,       .5+material.y1,     .5-material.z1;
                                                    -material.x1+1,       -material.y1+1,     -material.z1+1;
                                                       material.x2,          material.y2,        material.z2;
                                                    -material.x2+1,     .5+material.y2-1,     .5-material.z2;
                                                    -material.x2+1,       -material.y2+1,     -material.z2+1;
                                                    -material.x1+1,       -material.y1+1,     -material.z1+1;
                                                       material.x2,          material.y2,        material.z2;
                                                       material.x2,     .5-material.y2+1,     .5+material.z2;
                                                       material.x1,          material.y1,        material.z1;
                                                       material.x1,       .5-material.y1,     .5+material.z1;
                                                       material.x2,          material.y2,        material.z2;
                                                    -material.x2+1,       -material.y2+1,     -material.z2+1;
                                                       material.x2,     .5-material.y2+1,     .5+material.z2;
                                                    ];
material.parameters{1}.cylinder_top_centers = [          .5*(material.x1+material.x3),             .5*(material.y1+material.y3),           .5*(material.z1+material.z3);
                                                      .5*(-material.x2+1+material.x3),        .5*(.5+material.y2-1+material.y3),        .5*(.5-material.z2+material.z3);
                                                      .5*(-material.x2+1+material.x3),          .5*(-material.y2+1+material.y3),        .5*(-material.z2+1+material.z3);
                                                    .5*(-material.x1+1-material.x3+1),       .5*(.5+material.y1+.5+material.y3),     .5*(.5-material.z1+.5-material.z3);
                                                    .5*(-material.x1+1-material.x3+1),       .5*(-material.y1+1+.5+material.y3),     .5*(-material.z1+1+.5-material.z3);
                                                       .5*(material.x2-material.x3+1),          .5*(material.y2+.5+material.y3),        .5*(material.z2+.5-material.z3);
                                                    .5*(-material.x2+1-material.x3+1),     .5*(.5+material.y2-1+.5+material.y3),     .5*(.5-material.z2+.5-material.z3);
                                                    .5*(-material.x2+1-material.x3+1),       .5*(-material.y2+1+.5+material.y3),     .5*(-material.z2+1+.5-material.z3);
                                                    .5*(-material.x1+1-material.x3+1),        .5*(-material.y1+1-material.y3+1),      .5*(-material.z1+1-material.z3+1);
                                                       .5*(material.x2-material.x3+1),           .5*(material.y2-material.y3+1),         .5*(material.z2-material.z3+1);
                                                       .5*(material.x2-material.x3+1),      .5*(.5-material.y2+1-material.y3+1),      .5*(.5+material.z2-material.z3+1);
                                                         .5*(material.x1+material.x3),          .5*(material.y1+.5-material.y3),        .5*(material.z1+.5+material.z3);
                                                         .5*(material.x1+material.x3),       .5*(.5-material.y1+.5-material.y3),     .5*(.5+material.z1+.5+material.z3);
                                                         .5*(material.x2+material.x3),          .5*(material.y2+.5-material.y3),        .5*(material.z2+.5+material.z3);
                                                      .5*(-material.x2+1+material.x3),       .5*(-material.y2+1+.5-material.y3),     .5*(-material.z2+1+.5+material.z3);
                                                         .5*(material.x2+material.x3),     .5*(.5-material.y2+1+.5-material.y3),     .5*(.5+material.z2+.5+material.z3);
                                                    ];
material.parameters{1}.cylinder_bot_centers = material.parameters{1}.cylinder_bot_centers(:,[2,1,3]);
material.parameters{1}.cylinder_top_centers = material.parameters{1}.cylinder_top_centers(:,[2,1,3]);
material.parameters{1}.cylinder_radius      = cylinder_radius(1)*ones(size(material.parameters{1}.cylinder_top_centers,1),1);

%% Shape description for Zr  (Cylinders)
material.parameters{2}.cylinder_bot_centers = [        material.x3,        material.y3,        material.z3;
                                                       material.x3,        material.y3,        material.z3;
                                                       material.x3,        material.y3,        material.z3;
                                                    -material.x3+1,     .5+material.y3,     .5-material.z3;
                                                    -material.x3+1,     .5+material.y3,     .5-material.z3;
                                                    -material.x3+1,     .5+material.y3,     .5-material.z3;
                                                    -material.x3+1,     .5+material.y3,     .5-material.z3;
                                                    -material.x3+1,     .5+material.y3,     .5-material.z3;
                                                    -material.x3+1,     -material.y3+1,     -material.z3+1;
                                                    -material.x3+1,     -material.y3+1,     -material.z3+1;
                                                    -material.x3+1,     -material.y3+1,     -material.z3+1;
                                                       material.x3,     .5-material.y3,     .5+material.z3;
                                                       material.x3,     .5-material.y3,     .5+material.z3;
                                                       material.x3,     .5-material.y3,     .5+material.z3;
                                                       material.x3,     .5-material.y3,     .5+material.z3;
                                                       material.x3,     .5-material.y3,     .5+material.z3;
                                                    ];
material.parameters{2}.cylinder_top_centers = [          .5*(material.x1+material.x3),             .5*(material.y1+material.y3),           .5*(material.z1+material.z3);
                                                      .5*(-material.x2+1+material.x3),        .5*(.5+material.y2-1+material.y3),        .5*(.5-material.z2+material.z3);
                                                      .5*(-material.x2+1+material.x3),          .5*(-material.y2+1+material.y3),        .5*(-material.z2+1+material.z3);
                                                    .5*(-material.x1+1-material.x3+1),       .5*(.5+material.y1+.5+material.y3),     .5*(.5-material.z1+.5-material.z3);
                                                    .5*(-material.x1+1-material.x3+1),       .5*(-material.y1+1+.5+material.y3),     .5*(-material.z1+1+.5-material.z3);
                                                       .5*(material.x2-material.x3+1),          .5*(material.y2+.5+material.y3),        .5*(material.z2+.5-material.z3);
                                                    .5*(-material.x2+1-material.x3+1),     .5*(.5+material.y2-1+.5+material.y3),     .5*(.5-material.z2+.5-material.z3);
                                                    .5*(-material.x2+1-material.x3+1),       .5*(-material.y2+1+.5+material.y3),     .5*(-material.z2+1+.5-material.z3);
                                                    .5*(-material.x1+1-material.x3+1),        .5*(-material.y1+1-material.y3+1),      .5*(-material.z1+1-material.z3+1);
                                                       .5*(material.x2-material.x3+1),           .5*(material.y2-material.y3+1),         .5*(material.z2-material.z3+1);
                                                       .5*(material.x2-material.x3+1),      .5*(.5-material.y2+1-material.y3+1),      .5*(.5+material.z2-material.z3+1);
                                                         .5*(material.x1+material.x3),          .5*(material.y1+.5-material.y3),        .5*(material.z1+.5+material.z3);
                                                         .5*(material.x1+material.x3),       .5*(.5-material.y1+.5-material.y3),     .5*(.5+material.z1+.5+material.z3);
                                                         .5*(material.x2+material.x3),          .5*(material.y2+.5-material.y3),        .5*(material.z2+.5+material.z3);
                                                      .5*(-material.x2+1+material.x3),       .5*(-material.y2+1+.5-material.y3),     .5*(-material.z2+1+.5+material.z3);
                                                         .5*(material.x2+material.x3),     .5*(.5-material.y2+1+.5-material.y3),     .5*(.5+material.z2+.5+material.z3);
                                                    ];
material.parameters{2}.cylinder_bot_centers = material.parameters{2}.cylinder_bot_centers(:,[2,1,3]);
material.parameters{2}.cylinder_top_centers = material.parameters{2}.cylinder_top_centers(:,[2,1,3]);
material.parameters{2}.cylinder_radius      = cylinder_radius(2)*ones(size(material.parameters{2}.cylinder_top_centers,1),1);

%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};