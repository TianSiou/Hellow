function material = No227_FCC_NaTl(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No227_FCC_NaTl
%
% This is an example of a Zintl Phase.
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/AB_cF16_227_a_b.html
%
% Edit at 2017/7/17 By Hsiao-Han Huang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'NaTl';
material.material_num = 2;

material.lattice_type = 'face_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB_cF16_227_a_b.html';

material.lattice_constant.a = 7.483;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Na
material.parameters{1}.name                 = 'Na';
material.parameters{1}.color_map            = [128 0 255]/255;
material.parameters{1}.sphere_centers       = [     .125,     .125,     .125;
                                                    .875,     .875,     .875;
                                              ];                                               
material.parameters{1}.sphere_radius        =sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 

%% Shape description for Tl
material.parameters{2}.name                 = 'Tl';
material.parameters{2}.color_map            = [81 0 0]/255;
material.parameters{2}.sphere_centers       = [     .375,     .375,     .375;
                                                    .625,     .625,     .625;
                                              ];
                                         
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 

%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};