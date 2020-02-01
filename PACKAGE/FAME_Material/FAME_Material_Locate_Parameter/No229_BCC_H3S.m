function material = No229_BCC_H3S(sphere_radius, cylinder_radius)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No229_BCC_H3S
%
% See detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A3B_cI8_229_b_a.html
% Edit at 2017/6/28 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
material.Prototype    = 'H3S';
material.material_num = 2;

material.lattice_type = 'body_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A3B_cI8_229_b_a.html';

material.lattice_constant.a = 2.98400;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for S
material.parameters{1}.name                 = 'S';
material.parameters{1}.color_map            = [255 242 0]/255;
material.parameters{1}.sphere_centers       = [    0 ,   0 ,    0 ];
                                                  
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 
material.parameters{1}.cylinder_bot_centers =  [    0 ,   0 ,    0 ;
                                                    0 ,   0 ,    0 ;
                                                    0 ,   0 ,    0 ;
                                                    0 ,   1 ,    1 ;
                                                    1 ,   1 ,    1 ;
                                                    1 ,   1 ,    1 ];
                                               
material.parameters{1}.cylinder_top_centers = 0.5* [   0 ,    0.5 ,    0.5 ;  
                                                     0.5 ,      0 ,    0.5 ;  
                                                     0.5 ,    0.5 ,      0 ;
                                                       0 ,    1.5 ,    1.5 ;
                                                     1.5 ,    1.5 ,      2 ;
                                                     1.5 ,      2 ,    1.5 ;]; 
material.parameters{1}.cylinder_radius        = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_bot_centers,1)); 
size(material.parameters{1}.cylinder_bot_centers )
size(material.parameters{1}.cylinder_top_centers)
size(material.parameters{1}.cylinder_radius)
%% Shape description for H
material.parameters{2}.name                 = 'H';
material.parameters{2}.color_map            = [255 255 255]/255;
material.parameters{2}.sphere_centers       = [     0 ,      0.5 ,        0.5 ;  
                                                  0.5 ,        0 ,        0.5 ;  
                                                  0.5 ,      0.5 ,          0 ];                                    
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 
material.parameters{2}.cylinder_bot_centers =  [material.parameters{2}.sphere_centers;
                                                   0 , 0.5 , 0.5;  
                                                   0.5 , 0.5 , 1;
                                                   0.5 , 1 , 0.5];
                                                                                                   
material.parameters{2}.cylinder_top_centers =  material.parameters{1}.cylinder_top_centers ;
material.parameters{2}.cylinder_radius = cylinder_radius(2)*ones(1,size(material.parameters{2}.cylinder_bot_centers,1)); 

%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};