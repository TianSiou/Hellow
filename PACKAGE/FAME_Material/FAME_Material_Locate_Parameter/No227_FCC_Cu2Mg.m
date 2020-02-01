function material = No227_FCC_Cu2Mg(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No227_FCC_Cu2Mg
%
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A2B_cF24_227_d_a.html
%
% Edit at 2017/7/17 By Hsiao-Han Huang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'Cu2Mg';
material.material_num = 2;

material.lattice_type = 'face_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A2B_cF24_227_d_a.html';

material.lattice_constant.a = 7.02;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Mg
material.parameters{1}.color_map            = [0 255 0]/255;
material.parameters{1}.name                 = 'Mg';
material.parameters{1}.sphere_centers       = [     .125,     .125,     .125;
                                                    .875,     .875,     .875;
                                              ];                                         
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 

%% Shape description for Cu
material.parameters{2}.name                 = 'Cu';
material.parameters{2}.color_map            = [145 72 0]/255;
material.parameters{2}.sphere_centers       = [     .5,     .5,     .5;
                                                    .5,     .5,      0;
                                                    .5,      0,     .5;
                                                     0,     .5,     .5;
                                              ];
radius = 0.06*material.lattice_constant.a;                                               
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 

%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};