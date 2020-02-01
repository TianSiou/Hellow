function material = No224_SC_Cu2O(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No224_SC_Cu2O
%
% (Restori, 1986) gives the equilibrium lattice constant of Cu2O as a=4.627\AA, but gives nearest-neighbor distances which yield a lattice constant of 4.267\AA. 
% Since this value agrees with other sources, including those in (Downs, 2003), we use it.
%        http://www.aflowlib.org/CrystalDatabase/A2B_cP6_224_b_a.html
%
% Edit at 2017/7/17 By Hsiao-Han Huang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'Cu2O';
material.material_num = 2;

material.lattice_type = 'simple_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A2B_cP6_224_b_a.html';

material.lattice_constant.a = 4.267;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for O
material.parameters{1}.name                 = 'O';
material.parameters{1}.color_map            = [255 0 0]/255;
material.parameters{1}.sphere_centers       = [     .25,     .25,     .25;
                                                    .75,     .75,     .75;
                                              ];
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 

material.parameters{1}.cylinder_bot_centers = [     .25,     .25,     .25;
                                                    .25,     .25,     .25;
                                                    .25,     .25,     .25;
                                                    .25,     .25,     .25;
                                                    .75,     .75,     .75;
                                                    .75,     .75,     .75;
                                                    .75,     .75,     .75;
                                                    .75,     .75,     .75;
                                              ]; 

material.parameters{1}.cylinder_top_centers = [      .5*.25,      .5*.25,      .5*.25; 
                                                     .5*.75,      .5*.75,      .5*.25;
                                                     .5*.75,      .5*.25,      .5*.75;
                                                     .5*.25,      .5*.75,      .5*.75;
                                                    .5*1.75,     .5*1.75,     .5*1.75;
                                                    .5*1.25,     .5*1.25,     .5*1.75;
                                                    .5*1.25,     .5*1.75,     .5*1.25;
                                                    .5*1.75,     .5*1.25,     .5*1.25;
                                              ];
material.parameters{1}.cylinder_radius        = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_bot_centers,1)); 

%% Shape description for Cu
material.parameters{2}.name                 = 'Cu';
material.parameters{2}.color_map            = [202 101 0]/255;
material.parameters{2}.sphere_centers       = [      0,      0,      0;
                                                    .5,     .5,      0;
                                                    .5,      0,     .5;
                                                     0,     .5,     .5;
                                              ];
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 

material.parameters{2}.cylinder_bot_centers = [      0,      0,      0;
                                                    .5,     .5,      0;
                                                    .5,      0,     .5;
                                                     0,     .5,     .5;
                                                     1,      1,      1;
                                                    .5,     .5,      1;
                                                    .5,      1,     .5;
                                                     1,     .5,     .5;
                                              ];
                    
material.parameters{2}.cylinder_top_centers = [      .5*.25,      .5*.25,      .5*.25;
                                                     .5*.75,      .5*.75,      .5*.25;
                                                     .5*.75,      .5*.25,      .5*.75;
                                                     .5*.25,      .5*.75,      .5*.75;
                                                    .5*1.75,     .5*1.75,     .5*1.75;
                                                    .5*1.25,     .5*1.25,     .5*1.75;
                                                    .5*1.25,     .5*1.75,     .5*1.25;
                                                    .5*1.75,     .5*1.25,     .5*1.25;
                                              ];
material.parameters{2}.cylinder_radius        = cylinder_radius(2)*ones(1,size(material.parameters{2}.cylinder_bot_centers,1));
%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};