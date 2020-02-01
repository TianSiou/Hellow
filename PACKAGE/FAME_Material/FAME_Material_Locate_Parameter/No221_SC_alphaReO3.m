function material = No221_SC_alphaReO3(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No221_SC_alphaReO3
%
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A3B_cP4_221_d_a.html
%
% Edit at 2017/7/11 By Hsiao-Han Huang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'alphaReO3';
material.material_num = 2;

material.lattice_type = 'simple_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A3B_cP4_221_d_a.html';

material.lattice_constant.a = 3.734;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Re
material.parameters{1}.name                 = 'Re';
material.parameters{1}.color_map            = [0 0 255]/255;
material.parameters{1}.sphere_centers       = [   0,   0,   0];                                              
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 

material.parameters{1}.cylinder_bot_centers = [   0,   0,   0;
                                                  0,   0,   0;
                                                  0,   0,   0;
                                                  0,   1,   0;
                                                  0,   0,   1;
                                                  1,   0,   0]; 
material.parameters{1}.cylinder_top_centers = [    .5*.5,        0,        0;
                                                       0,    .5*.5,        0;
                                                       0,        0,    .5*.5;
                                                       0,   .5*1.5,        0;
                                                       0,        0,   .5*1.5;
                                                  .5*1.5,        0,        0];
                                              
material.parameters{1}.cylinder_radius        = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_bot_centers,1)); 
%% Shape description for O
material.parameters{2}.name                 = 'O';
material.parameters{2}.color_map            = [255 0 0]/255;
material.parameters{2}.sphere_centers       = [   .5,    0,    0;
                                                   0,   .5,    0;
                                                   0,    0,   .5];
                                               
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 

material.parameters{2}.cylinder_bot_centers = [   .5,    0,    0;
                                                   0,   .5,    0;
                                                   0,    0,   .5;
                                                   0,   .5,    0;
                                                   0,    0,   .5;
                                                  .5,    0,    0;
                                                      ];
                                              
material.parameters{2}.cylinder_top_centers = [    .5*.5,        0,        0;
                                                       0,    .5*.5,        0;
                                                       0,        0,    .5*.5;
                                                       0,   .5*1.5,        0;
                                                       0,        0,   .5*1.5;
                                                  .5*1.5,        0,        0;
                                                   ];
                                           
material.parameters{2}.cylinder_radius        = cylinder_radius(2)*ones(1,size(material.parameters{2}.cylinder_bot_centers,1)); 
%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};
