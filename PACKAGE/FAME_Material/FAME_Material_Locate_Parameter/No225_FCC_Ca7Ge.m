function material = No225_FCC_Ca7Ge(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No225_FCC_Ca7Ge
%
% See detail on following website:
%      http://www.aflowlib.org/CrystalDatabase/A7B_cF32_225_bd_a.html
% Edit at 2017/6/20 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'Ca7Ge';
material.material_num = 2;

material.lattice_type = 'face_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A7B_cF32_225_bd_a.html';
material.lattice_constant.a = 9.45000;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Ge
material.parameters{1}.name                 = 'Ge';
material.parameters{1}.color_map            = [0 128 64]/255;
material.parameters{1}.sphere_centers       = [    0 , 0 , 0 ] ;
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1));
    
%% Shape description for Ca
material.parameters{2}.name                 = 'Ca';
material.parameters{2}.color_map            = [0 255 0]/255;
material.parameters{2}.sphere_centers       = [    0.5 , 0.5 , 0.5 ;
                                                   0.5 ,   0 ,   0 ;
                                                     0 , 0.5 , 0.5 ;
                                                     0 , 0.5 ,   0 ;
                                                   0.5 ,   0 , 0.5 ;
                                                     0 ,   0 , 0.5 ;
                                                   0.5 , 0.5 ,   0 ];
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1));
    


%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};

    