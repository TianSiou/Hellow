function material = No225_FCC_NaCl(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No225_FCC_NaCl
%
% See detail on following website:
%      http://www.aflowlib.org/CrystalDatabase/AB_cF8_225_a_b.html
% Edit at 2017/6/22 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'NaCl';
material.material_num = 2;

material.lattice_type = 'face_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB_cF8_225_a_b.html';
material.lattice_constant.a = 5.63931;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Na
material.parameters{1}.name                 = 'Na';
material.parameters{1}.color_map            = [128 0 255]/255;
material.parameters{1}.sphere_centers       = [    0 ,  0 ,   0 ];                                          
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1));  
   
%% Shape description for Cl
material.parameters{2}.name                 = 'Cl';
material.parameters{2}.color_map            =  [0 255 0]/255;
material.parameters{2}.sphere_centers       = [  0.5 , 0.5 , 0.5  ];                                         
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1));  
   
%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map };