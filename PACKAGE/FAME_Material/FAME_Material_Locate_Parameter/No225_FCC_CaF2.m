function material = No225_FCC_CaF2(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No225_FCC_CaF2
%
% See detail on following website:
%      http://www.aflowlib.org/CrystalDatabase/AB2_cF12_225_a_c.html
% Edit at 2017/6/22 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'CaF2';
material.material_num = 2;

material.lattice_type = 'face_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB2_cF12_225_a_c.html';
material.lattice_constant.a = 5.46310;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Ca
material.parameters{1}.name                 = 'Ca';
material.parameters{1}.color_map            = [0 255 0]/255;
material.parameters{1}.sphere_centers       = [    0 ,  0 ,  0 ];
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1));
  

%% Shape description for F
material.parameters{2}.name                 = 'F';
material.parameters{2}.color_map            = [99 126 14]/255;
material.parameters{2}.sphere_centers       = [    0.25 ,  0.25 ,  0.25 ;
                                                   0.75 ,  0.75 ,  0.75 ];
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1));



%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};

    