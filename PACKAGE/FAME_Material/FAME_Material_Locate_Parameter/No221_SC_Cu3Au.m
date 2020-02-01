function material = No221_SC_Cu3Au(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No221_SC_Cu3Au
%
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/AB3_cP4_221_a_c.html
%
% Edit at 2017/7/17 By Hsiao-Han Huang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
material.Prototype    = 'Cu3Au';
material.material_num = 2;

material.lattice_type = 'simple_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB3_cP4_221_a_c.html';

material.lattice_constant.a = 3.7402;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Au
material.parameters{1}.name                 = 'Au';
material.parameters{1}.color_map            = [221 221 0]/255;
material.parameters{1}.sphere_centers       = [     0,     0,     0;
                                              ];
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 


%% Shape description for Cu
material.parameters{2}.name                 = 'Cu';
material.parameters{2}.color_map            = [153 77 0]/255;
material.parameters{2}.sphere_centers       = [      0,     .5,     .5;
                                                    .5,      0,     .5;
                                                    .5,     .5,      0;
                                              ];
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 

 
%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};