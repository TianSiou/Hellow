function material = No223_SC_Cr3Si(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No223_SC_Cr3Si
%
% The A Strukturbericht designation comes from the fact that this is also the structure of £]£]¡VW.
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A3B_cP8_223_c_a.html
%
% Edit at 2017/7/17 By Hsiao-Han Huang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'Cr3Si';
material.material_num = 2;

material.lattice_type = 'simple_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A3B_cP8_223_c_a.html';

material.lattice_constant.a = 4.556;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Si
material.parameters{1}.name                 = 'Si';
material.parameters{1}.color_map            = [255 179 102]/255;
material.parameters{1}.sphere_centers       = [      0,      0,      0;
                                                    .5,     .5,     .5;
                                              ];
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 

material.parameters{1}.cylinder_bot_centers = [      0,      0,      0;
                                                     0,      0,      0;
                                                     0,      0,      0;
                                                     1,      1,      0;
                                                     1,      1,      0;
                                                     1,      1,      0;
                                                     1,      0,      1;
                                                     1,      0,      1;
                                                     1,      0,      1;
                                                     0,      1,      1;
                                                     0,      1,      1;
                                                     0,      1,      1;
                                                    .5,     .5,     .5;
                                                    .5,     .5,     .5;
                                                    .5,     .5,     .5;
                                                    .5,     .5,     .5;
                                                    .5,     .5,     .5;
                                                    .5,     .5,     .5;
                                                    .5,     .5,     .5;
                                                    .5,     .5,     .5;
                                                    .5,     .5,     .5;
                                                    .5,     .5,     .5;
                                                    .5,     .5,     .5;
                                                    .5,     .5,     .5;
                                              ]; 

material.parameters{1}.cylinder_top_centers = [      .5*.25,           0,       .5*.5; 
                                                      .5*.5,      .5*.25,           0;
                                                          0,       .5*.5,      .5*.25;
                                                       .5*2,      .5*1.5,      .5*.25;
                                                     .5*1.5,     .5*1.75,           0;
                                                    .5*1.75,        .5*2,       .5*.5;
                                                     .5*1.5,      .5*.25,        .5*2;
                                                    .5*1.75,           0,      .5*1.5;
                                                       .5*2,       .5*.5,     .5*1.75;
                                                          0,      .5*1.5,     .5*1.75;
                                                      .5*.5,     .5*1.75,        .5*2;
                                                     .5*.25,        .5*2,      .5*1.5;
                                                     .5*.75,       .5*.5,          .5;
                                                    .5*1.25,       .5*.5,          .5;
                                                         .5,      .5*.75,       .5*.5;
                                                         .5,     .5*1.25,       .5*.5;
                                                      .5*.5,          .5,      .5*.75;
                                                      .5*.5,          .5,     .5*1.25;
                                                     .5*1.5,          .5,      .5*.75;
                                                    .5*1.25,      .5*1.5,          .5;
                                                         .5,      .5*.75,      .5*1.5;
                                                     .5*1.5,          .5,     .5*1.25;
                                                         .5,     .5*1.25,      .5*1.5;
                                                     .5*.75,      .5*1.5,          .5;
                                              ];
material.parameters{1}.cylinder_radius        = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_bot_centers,1)); 
%% Shape description for Cr
material.parameters{2}.name                 = 'Cr';
material.parameters{2}.color_map            = [132 132 255]/255;
material.parameters{2}.sphere_centers       = [     .25,       0,      .5;
                                                    .75,       0,      .5;
                                                     .5,     .25,       0;
                                                     .5,     .75,       0;
                                                      0,      .5,     .25;
                                                      0,      .5,     .75;
                                              ];
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 

material.parameters{2}.cylinder_bot_centers = [      .25,       0,      .5;
                                                      .5,     .25,       0;
                                                       0,      .5,     .25;
                                                       1,      .5,     .25;
                                                      .5,     .75,       0;
                                                     .75,       1,      .5;
                                                      .5,     .25,       1;
                                                     .75,       0,      .5;
                                                       1,      .5,     .75;
                                                       0,      .5,     .75;
                                                      .5,     .75,       1;
                                                     .25,       1,      .5;
                                                     .25,       0,      .5;
                                                     .75,       0,      .5;
                                                      .5,     .25,       0;
                                                      .5,     .75,       0;
                                                       0,      .5,     .25;
                                                       0,      .5,     .75;
                                                       1,      .5,     .25;
                                                     .75,       1,      .5;
                                                      .5,     .25,       1;
                                                       1,      .5,     .75;
                                                      .5,     .75,       1;
                                                     .25,       1,      .5;
                                              ];
                    
material.parameters{2}.cylinder_top_centers = [       .5*.25,           0,       .5*.5; 
                                                       .5*.5,      .5*.25,           0;
                                                           0,       .5*.5,      .5*.25;
                                                        .5*2,      .5*1.5,      .5*.25;
                                                      .5*1.5,     .5*1.75,           0;
                                                     .5*1.75,        .5*2,       .5*.5;
                                                      .5*1.5,      .5*.25,        .5*2;
                                                     .5*1.75,           0,      .5*1.5;
                                                        .5*2,       .5*.5,     .5*1.75;
                                                           0,      .5*1.5,     .5*1.75;
                                                       .5*.5,     .5*1.75,        .5*2;
                                                      .5*.25,        .5*2,      .5*1.5;
                                                      .5*.75,       .5*.5,          .5;
                                                     .5*1.25,       .5*.5,          .5;
                                                          .5,      .5*.75,       .5*.5;
                                                          .5,     .5*1.25,       .5*.5;
                                                       .5*.5,          .5,      .5*.75;
                                                       .5*.5,          .5,     .5*1.25;
                                                      .5*1.5,          .5,      .5*.75;
                                                     .5*1.25,      .5*1.5,          .5;
                                                          .5,      .5*.75,      .5*1.5;
                                                      .5*1.5,          .5,     .5*1.25;
                                                          .5,     .5*1.25,      .5*1.5;
                                                      .5*.75,      .5*1.5,          .5;
                                              ];
material.parameters{2}.cylinder_radius        = cylinder_radius(2)*ones(1,size(material.parameters{2}.cylinder_bot_centers,1)); 
%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};