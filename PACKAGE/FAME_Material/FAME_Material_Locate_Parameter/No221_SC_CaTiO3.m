function material = No221_SC_CaTiO3(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No221_SC_CaTiO3
%
% Cubic perovskite is actually the high-temperature phase of the compounds listed below. 
% The ground states are usually distorted perovskite structures. Many of these substances are ferroelectric. 
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/AB3C_cP5_221_a_c_b.html
%
% Edit at 2017/7/7 By Hsiao-Han Huang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'CaTiO3';
material.material_num = 3;

material.lattice_type = 'simple_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB3C_cP5_221_a_c_b.html';

material.lattice_constant.a = 3.795;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Ca
material.parameters{1}.name                 = 'Ca';
material.parameters{1}.color_map            = [46 141 57]/255;
material.parameters{1}.sphere_centers       = [   0,   0,   0];
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 

%% Shape description for Ti
material.parameters{2}.name                 = 'Ti';
material.parameters{2}.color_map            = [68 187 213]/255;
material.parameters{2}.sphere_centers       = [   .5,   .5,   .5];
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 

material.parameters{2}.cylinder_bot_centers = [   .5,   .5,   .5;
                                                  .5,   .5,   .5;
                                                  .5,   .5,   .5;
                                                  .5,   .5,   .5;
                                                  .5,   .5,   .5;
                                                  .5,   .5,   .5;
                                                      ];
                                              
material.parameters{2}.cylinder_top_centers = [   .5*.5,       .5*1,        .5*1;
                                                   .5*1,      .5*.5,        .5*1;
                                                   .5*1,       .5*1,       .5*.5;
                                               .5*(1.5),       .5*1,        .5*1;
                                                   .5*1,   .5*(1.5),        .5*1;
                                                   .5*1,       .5*1,    .5*(1.5)];
material.parameters{2}.cylinder_radius        = cylinder_radius(2)*ones(1,size(material.parameters{2}.cylinder_bot_centers,1)); 
%% Shape description for O
material.parameters{3}.name                 = 'O';
material.parameters{3}.color_map            = [233 13 1]/255;
material.parameters{3}.sphere_centers       = [   0,   .5,   .5;
                                                 .5,    0,   .5;
                                                 .5,   .5,    0];
material.parameters{3}.sphere_radius        = sphere_radius(3)*ones(1,size(material.parameters{3}.sphere_centers,1)); 

material.parameters{3}.cylinder_bot_centers = [   0,   .5,   .5;
                                                 .5,    0,   .5;
                                                 .5,   .5,    0;
                                                  1,   .5,   .5;
                                                 .5,    1,   .5;
                                                 .5,   .5,    1];
                                             
material.parameters{3}.cylinder_top_centers = [  .5*.5,       .5*1,        .5*1;
                                                  .5*1,      .5*.5,        .5*1;
                                                  .5*1,       .5*1,       .5*.5;
                                              .5*(1.5),       .5*1,        .5*1;
                                                  .5*1,   .5*(1.5),        .5*1;
                                                  .5*1,       .5*1,    .5*(1.5)];
material.parameters{3}.cylinder_radius        = cylinder_radius(3)*ones(1,size(material.parameters{3}.cylinder_bot_centers,1)); 
%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map, material.parameters{3}.color_map};