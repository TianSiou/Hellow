function material = No229_BCC_W(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No229_BCC_W
%
% Although more accurate measurements of the lattice constant of tungsten are available, (Davey, 1925) is chosen because of the unique experimental technique.
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A_cI2_229_a.html
%
% Edit at 2017/7/17 By Hsiao-Han Huang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'W';
material.material_num = 1;

material.lattice_type = 'body_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A_cI2_229_a.html';

material.lattice_constant.a = 3.155;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for W
material.parameters{1}.name                 = 'W';
material.parameters{1}.color_map            = [0 130 191]/255;

material.parameters{1}.sphere_centers       = [     0,     0,     0];
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 

%% Collect color map
material.color_map = {material.parameters{1}.color_map};