function material = No221_SC_CsCl(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No221_SC_CsCl
%
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/AB_cP2_221_b_a.html
%
% Edit at 2017/7/17 By Hsiao-Han Huang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
material.Prototype    = 'CsCl';
material.material_num = 2;

material.lattice_type = 'simple_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB_cP2_221_b_a.html';

material.lattice_constant.a = 4.07925;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Cs
material.parameters{1}.name                 = 'Cs';
material.parameters{1}.color_map            = [64 0 128]/255;
material.parameters{1}.sphere_centers       = [     0,     0,     0;
                                              ];
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 


%% Shape description for Cl
material.parameters{2}.name                 = 'Cl';
material.parameters{2}.color_map            = [0 255 0]/255;
material.parameters{2}.sphere_centers       = [     .5,     .5,     .5;
                                              ];
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 

 
%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};