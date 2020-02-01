function material = No221_SC_CaB6(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No221_SC_CaB6
%
% See detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A6B_cP7_221_f_a.html
% Edit at 2017/7/7 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
material.Prototype    = 'CaB6';
material.material_num = 2;

material.lattice_type = 'simple_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A6B_cP7_221_f_a.html';
material.x2           = 0.2117 ;
material.lattice_constant.a = 4.14500;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Ca
material.parameters{1}.name                 = 'Ca';
material.parameters{1}.color_map            = [0 255 0]/255;
material.parameters{1}.sphere_centers       = [    0 ,   0 ,   0 ];
                                      
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1));  
%% Shape description for B
material.parameters{2}.name                 = 'B';
material.parameters{2}.color_map            = [255 174 201]/255;
material.parameters{2}.sphere_centers       = [  material.x2  ,            0.5 ,           0.5 ; %2
                                                -material.x2  ,            0.5 ,           0.5 ; %3
                                                          0.5 ,   material.x2  ,           0.5 ; %4
                                                          0.5 ,  -material.x2  ,           0.5 ; %5
                                                          0.5 ,            0.5 ,   material.x2 ; %6
                                                          0.5 ,            0.5 ,  -material.x2 ]; %7
                                                                                           
material.parameters{2}.sphere_centers = mod(material.parameters{2}.sphere_centers,1);
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 
material.parameters{2}.cylinder_bot_centers =  [            0.5 ,            0.5 ,   material.x2 ; %6
                                                            0.5 ,            0.5 ,   material.x2 ; %6
                                                            0.5 ,            0.5 ,   material.x2 ; %6
                                                            0.5 ,            0.5 ,   material.x2 ; %6
                                                            0.5 ,            0.5 ,  -material.x2 ; %7
                                                            0.5 ,            0.5 ,  -material.x2 ; %7
                                                            0.5 ,            0.5 ,  -material.x2 ; %7
                                                            0.5 ,            0.5 ,  -material.x2 ; %7
                                                   material.x2  ,            0.5 ,           0.5 ; %2
                                                   material.x2  ,            0.5 ,           0.5 ; %2
                                                  -material.x2  ,            0.5 ,           0.5 ; %3
                                                  -material.x2  ,            0.5 ,           0.5 ]; %3
                                                     
material.parameters{2}.cylinder_bot_centers = mod(material.parameters{2}.cylinder_bot_centers,1);                                                                                                   
material.parameters{2}.cylinder_top_centers =  [   material.x2  ,            0.5 ,           0.5 ; %2
                                                  -material.x2  ,            0.5 ,           0.5 ; %3
                                                            0.5 ,   material.x2  ,           0.5 ; %4
                                                            0.5 ,  -material.x2  ,           0.5 ; %5
                                                   material.x2  ,            0.5 ,           0.5 ; %2
                                                  -material.x2  ,            0.5 ,           0.5 ; %3
                                                            0.5 ,   material.x2  ,           0.5 ; %4
                                                            0.5 ,  -material.x2  ,           0.5 ; %5
                                                            0.5 ,   material.x2  ,           0.5 ; %4
                                                            0.5 ,  -material.x2  ,           0.5 ; %5
                                                            0.5 ,   material.x2  ,           0.5 ; %4
                                                            0.5 ,  -material.x2  ,           0.5 ]; %5
material.parameters{2}.cylinder_top_centers = mod(material.parameters{2}.cylinder_top_centers,1);       
material.parameters{2}.cylinder_radius = cylinder_radius(2)*ones(1,size(material.parameters{2}.cylinder_bot_centers,1)); 


%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};