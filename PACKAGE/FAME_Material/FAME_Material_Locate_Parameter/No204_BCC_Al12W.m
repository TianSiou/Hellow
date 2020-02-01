function material = No204_BCC_Al12W(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No204_BCC_Al12W
%
% See detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A12B_cI26_204_g_a.html
%
% Edit at 2017/6/20 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'Al12W';
material.material_num = 2;

material.lattice_type = 'body_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A12B_cI26_204_g_a.html';
material.y2           = 0.184 ;
material.z2           = 0.691 ;

material.lattice_constant.a = 7.58000;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Ga
material.parameters{1}.name                 = 'W';
material.parameters{1}.color_map            = [0 128 255]/255;
material.parameters{1}.sphere_centers       = [       0 ,     0 ,       0 ];                       
material.parameters{1}.sphere_radius        = sphere_radius(1)*0.07;

%% Shape description for H
material.parameters{2}.name                 = 'Al';
material.parameters{2}.color_map            = [128 128 128]/255;
material.parameters{2}.sphere_centers       =  [    material.y2+material.z2 ,               material.z2 ,               material.y2 ;
                                                    material.z2-material.y2 ,               material.z2 ,             1-material.y2 ;
                                                  material.y2-material.z2+1 ,             1-material.z2 ,               material.y2 ;
                                                  1-material.y2-material.z2 ,             1-material.z2 ,             1-material.y2 ; 
                                                                material.y2 ,   material.y2+material.z2 ,               material.z2 ; 
                                                              1-material.y2 ,   material.z2-material.y2 ,               material.z2 ; 
                                                                material.y2 , 1+material.y2-material.z2 ,             1-material.z2 ;
                                                              1-material.y2 , 1-material.y2-material.z2 ,             1-material.z2 ; 
                                                                material.z2 ,               material.y2 ,   material.y2+material.z2 ; 
                                                                material.z2 ,             1-material.y2 ,   material.z2-material.y2 ; 
                                                              1-material.z2 ,               material.y2 , material.y2-material.z2+1 ; 
                                                              1-material.z2 ,             1-material.y2 ,1-material.y2-material.z2] ;
                                             
material.parameters{2}.sphere_radius = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1));                                                             

%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};
