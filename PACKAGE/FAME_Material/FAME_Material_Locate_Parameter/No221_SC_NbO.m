function material = No221_SC_NbO(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No221_SC_NbO
%
% This is the NaCl (B1) structure with 25% ordered vacancies on both the Na and Cl sites.
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/AB_cP6_221_c_d.html
%
% Edit at 2017/7/11 By Hsiao-Han Huang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'NbO';
material.material_num = 2;

material.lattice_type = 'simple_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB_cP6_221_c_d.html';

material.lattice_constant.a = 4.2101;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Nb
material.parameters{1}.name                 = 'Nb';
material.parameters{1}.color_map            = [0 0 255]/255;
material.parameters{1}.sphere_centers       = [    0,   .5,   .5;
                                                  .5,    0,   .5;
                                                  .5,   .5,    0];                                            
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 

material.parameters{1}.cylinder_bot_centers = [    0,   .5,   .5;
                                                   0,   .5,   .5;
                                                   0,   .5,   .5;
                                                   0,   .5,   .5;
                                                  .5,    0,   .5;
                                                  .5,    0,   .5;
                                                  .5,    0,   .5;
                                                  .5,    0,   .5;
                                                  .5,   .5,    0;
                                                  .5,   .5,    0;
                                                  .5,   .5,    0;
                                                  .5,   .5,    0;
                                              ]; 

material.parameters{1}.cylinder_top_centers = [        0,     .5*1,    .5*.5;
                                                       0,    .5*.5,     .5*1; 
                                                       0,   .5*1.5,     .5*1;
                                                       0,     .5*1,   .5*1.5;
                                                    .5*1,        0,    .5*.5;
                                                   .5*.5,        0,     .5*1;
                                                  .5*1.5,        0,     .5*1;
                                                    .5*1,        0,   .5*1.5;
                                                    .5*1,    .5*.5,        0;
                                                   .5*.5,     .5*1,        0;
                                                    .5*1,   .5*1.5,        0;
                                                  .5*1.5,     .5*1,        0;
                                              ];
material.parameters{1}.cylinder_radius        = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_bot_centers,1)); 
%% Shape description for O
material.parameters{2}.name                 = 'O';
material.parameters{2}.color_map            = [255 0 0]/255;
material.parameters{2}.sphere_centers       = [   .5,    0,    0;
                                                   0,   .5,    0;
                                                   0,    0,   .5];
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 

material.parameters{2}.cylinder_bot_centers = [    0,   .5,    0;
                                                   0,    0,   .5;
                                                   0,    1,   .5;
                                                   0,   .5,    1;
                                                  .5,    0,    0;
                                                   0,    0,   .5;
                                                   1,    0,   .5;
                                                  .5,    0,    1;
                                                  .5,    0,    0;
                                                   0,   .5,    0;
                                                  .5,    1,    0;
                                                   1,   .5,    0;
                                              ];
                    
material.parameters{2}.cylinder_top_centers = [        0,     .5*1,    .5*.5;
                                                       0,    .5*.5,     .5*1;
                                                       0,   .5*1.5,     .5*1;
                                                       0,     .5*1,   .5*1.5;
                                                    .5*1,        0,    .5*.5;
                                                   .5*.5,        0,     .5*1;
                                                  .5*1.5,        0,     .5*1;
                                                    .5*1,        0,   .5*1.5;
                                                    .5*1,    .5*.5,        0;
                                                   .5*.5,     .5*1,        0;
                                                    .5*1,   .5*1.5,        0;
                                                  .5*1.5,     .5*1,        0;
                                              ];
material.parameters{2}.cylinder_radius        = cylinder_radius(2)*ones(1,size(material.parameters{2}.cylinder_bot_centers,1)); 
%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};
