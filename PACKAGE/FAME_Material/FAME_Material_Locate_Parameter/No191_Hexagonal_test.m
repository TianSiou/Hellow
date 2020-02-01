function material = No191_Hexagonal_test(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No191_Hexagonal_test
%
% More detail on following website:
%        https://www.nature.com/articles/nphoton.2016.253.pdf
%
% Edit at 2018/1/16 By Jia-Wei Lin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'test';
material.material_num = 1;

material.lattice_type = 'hexagonal';
material.website      = 'https://www.nature.com/articles/nphoton.2016.253.pdf';

% r = 1e+8;
r = 1;
material.a0 = 6.5*r;
material.az = 0.69*material.a0;
material.hr = (2.11*r)/material.az;

material.lattice_constant.a = material.a0;
material.lattice_constant.c = material.az;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for C
material.parameters{1}.name                 = 'C';
material.parameters{1}.color_map            = [238 94 6]/255;
material.parameters{1}.cylinder_top_centers       = [   .5, .5, .5-material.hr/2];
material.parameters{1}.cylinder_bot_centers       = [   .5, .5, .5+material.hr/2];    
% material.parameters{1}.cylinder_top_centers       = [   0, 0, 0;
%                                                         1, 0, 0;
%                                                         0, 1, 0;
%                                                         1, 1, 0;
%                                                         0, 0, 1;
%                                                         1, 0, 1;
%                                                         0, 1, 1;
%                                                         1, 1, 1];
% material.parameters{1}.cylinder_bot_centers       = [   0, 0, material.hr/2;
%                                                         1, 0, material.hr/2;
%                                                         0, 1, material.hr/2;
%                                                         1, 1, material.hr/2;  
%                                                         0, 0, 1-material.hr/2;
%                                                         1, 0, 1-material.hr/2;
%                                                         0, 1, 1-material.hr/2;
%                                                         1, 1, 1-material.hr/2];                                                    
                                          
material.parameters{1}.cylinder_radius        = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_bot_centers,1)); 

%% Collect color map
material.color_map = {material.parameters{1}.color_map};