function material = No221_SC_alpha_Po(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No221_SC_alpha_Po
%
% See detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A_cP1_221_a.html
%
% Edit at 2017/6/16 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'No221_SC_alpha_Po';
material.material_num = 1;

material.lattice_type = 'simple_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A_cP1_221_a.html';

material.lattice_constant.a = 3.34;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Fe
material.parameters{1}.name                 = 'Po';
material.parameters{1}.color_map            = [168 109 57]/255;
material.parameters{1}.sphere_centers       = [ 0, 0, 0 ];                                               
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 

%% Collect color map
material.color_map = {material.parameters{1}.color_map};
