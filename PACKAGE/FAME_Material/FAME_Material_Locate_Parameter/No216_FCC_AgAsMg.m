function material = No216_FCC_AgAsMg(sphere_radius, cylinder_radius)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No216_FCC_AgAsMg
%
% See detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/ABC_cF12_216_b_c_a.html
% Edit at 2017/7/3 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
material.Prototype    = 'AgAsMg';
material.material_num = 3;

material.lattice_type = 'face_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/ABC_cF12_216_b_c_a.html';

material.lattice_constant.a = 6.24000;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Mg
material.parameters{1}.name                 = 'Mg';
material.parameters{1}.color_map            = [128 255 0]/255;
material.parameters{1}.sphere_centers       = [    0 ,   0 ,    0 ];
                                                  
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 
material.parameters{1}.cylinder_bot_centers =  [material.parameters{1}.sphere_centers;
                                                 0 , 0 , 1 ;
                                                 1 , 0 , 0 ;
                                                 0 , 1 , 0  ];
                                               
material.parameters{1}.cylinder_top_centers = 0.5* [  0.25 , 0.25 , 0.25 ;
                                                      0.25 , 0.25 , 1.25 ;
                                                      1.25 , 0.25 , 0.25 ;
                                                      0.25 , 1.25 , 0.25]; 
material.parameters{1}.cylinder_radius        = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_bot_centers,1)); 
%% Shape description for Ag
material.parameters{2}.name                 = 'Ag';
material.parameters{2}.color_map            = [195 195 195]/255;
material.parameters{2}.sphere_centers       = [     0.5 ,  0.5 ,    0.5];                                    

material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 
material.parameters{2}.cylinder_bot_centers =  [material.parameters{2}.sphere_centers;
                                                0.5 , -0.5 , 0.5 ; 
                                               -0.5 ,  0.5 , 0.5 ] ;
                                                                                                   
material.parameters{2}.cylinder_top_centers =  0.5*[0.75 , 0.75 , 0.75 ;
                                                    0.75 ,-0.25 , 0.75 ;
                                                   -0.25 , 0.75 , 0.75 ] ;
%material.parameters{2}.cylinder_bot_centers = mod(material.parameters{2}.cylinder_bot_centers,1);
%material.parameters{2}.cylinder_top_centers = mod(material.parameters{2}.cylinder_top_centers,1);
material.parameters{2}.cylinder_radius = cylinder_radius(2)*ones(1,size(material.parameters{2}.cylinder_bot_centers,1)); 

%% Shape description for As
material.parameters{3}.name                 = 'As';
material.parameters{3}.color_map            = [163 73 164]/255;
material.parameters{3}.sphere_centers       = [    0.25 ,    0.25 ,     0.25 ];                                    

material.parameters{3}.sphere_radius        = sphere_radius(3)*ones(1,size(material.parameters{2}.sphere_centers,1)); 
material.parameters{3}.cylinder_bot_centers =  [  0.25 ,    0.25 ,     0.25 ;
                                                  0.25 ,    0.25 ,     0.25 ;
                                                  0.25 ,    0.25 ,     0.25 ;
                                                  0.25 ,    0.25 ,     0.25 ;
                                                  0.25 ,    0.25 ,     0.25 ;
                                                  0.25 ,    0.25 ,     0.25 ;
                                                  0.25 ,    0.25 ,     0.25];
                                                                                                   
material.parameters{3}.cylinder_top_centers =  0.5*[  0.25 , 0.25 , 0.25 ;
                                                      0.75 , 0.75 , 0.75 ;
                                                      0.25 , 0.25 , 1.25 ;
                                                      1.25 , 0.25 , 0.25 ;
                                                      0.25 , 1.25 , 0.25 ;
                                                      0.75 , -0.25 , 0.75;
                                                      -0.25 , 0.75 , 0.75] ;
%material.parameters{3}.cylinder_top_centers = mod(material.parameters{3}.cylinder_top_centers,1); 
material.parameters{3}.cylinder_radius = cylinder_radius(3)*ones(1,size(material.parameters{3}.cylinder_bot_centers,1)); 
%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map, material.parameters{3}.color_map};