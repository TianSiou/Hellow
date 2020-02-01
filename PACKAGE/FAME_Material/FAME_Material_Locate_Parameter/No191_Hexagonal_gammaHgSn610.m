function material = No191_Hexagonal_gammaHgSn610(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No191_Hexagonal_gammaHgSn610
%
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A_hP1_191_a.html
%
% Edit at 2017/7/25 by Jyun-Wei Lin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'gammaHgSn610';
material.material_num = 1;

material.lattice_type = 'hexagonal';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A_hP1_191_a.html';

material.lattice_constant.a = 3.2062;
material.lattice_constant.c = 0.931195808122*material.lattice_constant.a;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Ca
material.parameters{1}.name                 = 'Sn';
material.parameters{1}.color_map            = [88 88 88]/255;
material.parameters{1}.sphere_centers       = [   0,  0,  0
                                              ];
                                             
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 

%% Collect color map
material.color_map = {material.parameters{1}.color_map};