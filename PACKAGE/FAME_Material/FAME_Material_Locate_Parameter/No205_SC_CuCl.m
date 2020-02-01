function material = No205_SC_CuCl(sphere_radius, cylinder_radius)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No205_SC_CuCl
%
% See detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/AB_cP16_205_c_c.html
% Edit at 2017/7/5 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
material.Prototype    = 'CuCl';
material.material_num = 2;

material.lattice_type = 'simple_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB_cP16_205_c_c.html';
material.x1           = 0.1527 ;
material.x2           = 0.6297 ;
material.lattice_constant.a = 6.41620;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Cu
material.parameters{1}.name                 = 'Cu';
material.parameters{1}.color_map            = [158 109 50]/255;
material.parameters{1}.sphere_centers       = [      material.x2 ,     material.x2 ,   material.x2   ;  %9  
                                                 0.5-material.x2 ,    -material.x2 , 0.5+material.x2 ;  %10
                                                    -material.x2 , 0.5+material.x2 , 0.5-material.x2 ;  %11
                                                 0.5+material.x2 , 0.5-material.x2 ,    -material.x2 ;  %12
                                                    -material.x2 ,    -material.x2 ,    -material.x2 ;  %13
                                                 0.5+material.x2 ,     material.x2 , 0.5-material.x2 ;  %14
                                                     material.x2 , 0.5-material.x2 , 0.5+material.x2 ;  %15
                                                 0.5-material.x2 , 0.5+material.x2 ,      material.x2]; %16 
material.parameters{1}.sphere_centers = mod(material.parameters{1}.sphere_centers,1);                                                  

material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 
material.parameters{1}.cylinder_bot_centers =  [    material.x2 ,     material.x2 ,   material.x2   ;  %9  
                                                    material.x2 ,     material.x2 ,   material.x2   ;  %9  
                                                    material.x2 ,     material.x2 ,   material.x2   ;  %9  
                                                    material.x2 ,     material.x2 ,   material.x2   ;  %9  
                                                0.5-material.x2 ,    -material.x2 , 0.5+material.x2 ;  %10
                                                0.5-material.x2 ,    -material.x2 , 0.5+material.x2 ;  %10
                                                   -material.x2 , 0.5+material.x2 , 0.5-material.x2 ;  %11
                                                   -material.x2 , 0.5+material.x2 , 0.5-material.x2 ;  %11
                                                0.5+material.x2 , 0.5-material.x2 ,    -material.x2 ;  %12
                                                0.5+material.x2 , 0.5-material.x2 ,    -material.x2 ;  %12
                                                   -material.x2 ,    -material.x2 ,    -material.x2 ;  %13
                                                   -material.x2 ,    -material.x2 ,    -material.x2 ;  %13
                                                   -material.x2 ,    -material.x2 ,    -material.x2 ;  %13
                                                   -material.x2 ,    -material.x2 ,    -material.x2 ;  %13
                                                0.5+material.x2 ,     material.x2 , 0.5-material.x2 ;  %14
                                                0.5+material.x2 ,     material.x2 , 0.5-material.x2 ;  %14
                                                    material.x2 , 0.5-material.x2 , 0.5+material.x2 ;  %15
                                                    material.x2 , 0.5-material.x2 , 0.5+material.x2 ;  %15
                                                0.5-material.x2 , 0.5+material.x2 ,     material.x2 ;  %16 
                                                0.5-material.x2 , 0.5+material.x2 ,     material.x2 ]; %16  
material.parameters{1}.cylinder_bot_centers = mod(material.parameters{1}.cylinder_bot_centers,1);                                               
material.parameters{1}.cylinder_top_centers =  [  0.5-material.x1 ,    -material.x1 , 0.5+material.x1 ; %2
                                                     -material.x1 , 0.5+material.x1 , 0.5-material.x1 ; %3
                                                  0.5+material.x1 , 0.5-material.x1 ,    -material.x1 ; %4
                                                     -material.x1 ,    -material.x1 ,    -material.x1 ; %5
                                                     -material.x1 , 0.5+material.x1 , 0.5-material.x1 ; %3
                                                  0.5+material.x1 ,     material.x1 , 0.5-material.x1 ; %6
                                                  0.5+material.x1 , 0.5-material.x1 ,    -material.x1 ; %4
                                                      material.x1 , 0.5-material.x1 , 0.5+material.x1 ; %7
                                                  0.5-material.x1 ,    -material.x1 , 0.5+material.x1 ; %2
                                                  0.5-material.x1 , 0.5+material.x1 ,     material.x1 ; %8
                                                      material.x1 ,     material.x1 ,     material.x1 ; %1
                                                  0.5+material.x1 ,     material.x1 , 0.5-material.x1 ; %6
                                                      material.x1 , 0.5-material.x1 , 0.5+material.x1 ; %7
                                                  0.5-material.x1 , 0.5+material.x1 ,     material.x1 ; %8
                                                  0.5-material.x1 ,    -material.x1 , 0.5+material.x1 ; %2
                                                      material.x1 , 0.5-material.x1 , 0.5+material.x1 ; %7
                                                     -material.x1 , 0.5+material.x1 , 0.5-material.x1 ; %3
                                                  0.5-material.x1 , 0.5+material.x1 ,     material.x1 ; %8 
                                                  0.5+material.x1 , 0.5-material.x1 ,    -material.x1 ; %4
                                               0.5+material.x1 ,     material.x1 , 0.5-material.x1   ]; %6 
material.parameters{1}.cylinder_top_centers = mod(material.parameters{1}.cylinder_top_centers,1);
material.parameters{1}.cylinder_top_centers =0.5*(material.parameters{1}.cylinder_top_centers +material.parameters{1}.cylinder_bot_centers);

material.parameters{1}.cylinder_radius        = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_bot_centers,1)); 


%% Shape description for Cl
material.parameters{2}.name                 = 'Cl';
material.parameters{2}.color_map            = [0 255 0]/255;
material.parameters{2}.sphere_centers       = [       material.x1 ,     material.x1 ,     material.x1 ; %1
                                                  0.5-material.x1 ,    -material.x1 , 0.5+material.x1 ; %2
                                                     -material.x1 , 0.5+material.x1 , 0.5-material.x1 ; %3
                                                  0.5+material.x1 , 0.5-material.x1 ,    -material.x1 ; %4
                                                     -material.x1 ,    -material.x1 ,    -material.x1 ; %5
                                                  0.5+material.x1 ,     material.x1 , 0.5-material.x1 ; %6
                                                      material.x1 , 0.5-material.x1 , 0.5+material.x1 ; %7
                                                  0.5-material.x1 , 0.5+material.x1 ,     material.x1 ];%8                                 
material.parameters{2}.sphere_centers  = mod(material.parameters{2}.sphere_centers ,1);

material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 
material.parameters{2}.cylinder_bot_centers =  [      material.x1 ,     material.x1 ,     material.x1 ; %1
                                                  0.5-material.x1 ,    -material.x1 , 0.5+material.x1 ; %2
                                                  0.5-material.x1 ,    -material.x1 , 0.5+material.x1 ; %2
                                                  0.5-material.x1 ,    -material.x1 , 0.5+material.x1 ; %2
                                                     -material.x1 , 0.5+material.x1 , 0.5-material.x1 ; %3
                                                     -material.x1 , 0.5+material.x1 , 0.5-material.x1 ; %3
                                                     -material.x1 , 0.5+material.x1 , 0.5-material.x1 ; %3
                                                  0.5+material.x1 , 0.5-material.x1 ,    -material.x1 ; %4
                                                  0.5+material.x1 , 0.5-material.x1 ,    -material.x1 ; %4
                                                  0.5+material.x1 , 0.5-material.x1 ,    -material.x1 ; %4
                                                     -material.x1 ,    -material.x1 ,    -material.x1 ; %5
                                                  0.5+material.x1 ,     material.x1 , 0.5-material.x1 ; %6
                                                  0.5+material.x1 ,     material.x1 , 0.5-material.x1 ; %6
                                                  0.5+material.x1 ,     material.x1 , 0.5-material.x1 ; %6
                                                      material.x1 , 0.5-material.x1 , 0.5+material.x1 ; %7
                                                      material.x1 , 0.5-material.x1 , 0.5+material.x1 ; %7
                                                      material.x1 , 0.5-material.x1 , 0.5+material.x1 ; %7
                                                  0.5-material.x1 , 0.5+material.x1 ,     material.x1 ; %8
                                                  0.5-material.x1 , 0.5+material.x1 ,     material.x1 ; %8
                                                  0.5-material.x1 , 0.5+material.x1 ,     material.x1 ];%8
material.parameters{2}.cylinder_bot_centers = mod(material.parameters{2}.cylinder_bot_centers,1);                                                                                                   
material.parameters{2}.cylinder_top_centers =  [          -material.x2 ,    -material.x2 ,    -material.x2 ;  %13
                                                           material.x2 ,     material.x2 ,   material.x2   ;  %9 
                                                       0.5+material.x2 , 0.5-material.x2 ,    -material.x2 ;  %12
                                                       0.5+material.x2 ,     material.x2 , 0.5-material.x2 ;  %14
                                                           material.x2 ,     material.x2 ,   material.x2   ;  %9  
                                                       0.5-material.x2 ,    -material.x2 , 0.5+material.x2 ;  %10
                                                           material.x2 , 0.5-material.x2 , 0.5+material.x2 ;  %15
                                                           material.x2 ,     material.x2 ,     material.x2 ;  %9   
                                                          -material.x2 , 0.5+material.x2 , 0.5-material.x2 ;  %11
                                                       0.5-material.x2 , 0.5+material.x2 ,      material.x2;  %16
                                                           material.x2 ,     material.x2 ,     material.x2 ;  %9 
                                                       0.5-material.x2 ,    -material.x2 , 0.5+material.x2 ;  %10
                                                          -material.x2 ,    -material.x2 ,    -material.x2 ;  %13
                                                       0.5-material.x2 , 0.5+material.x2 ,      material.x2;  %16
                                                          -material.x2 , 0.5+material.x2 , 0.5-material.x2 ;  %11 
                                                          -material.x2 ,    -material.x2 ,    -material.x2 ;  %13
                                                       0.5+material.x2 ,     material.x2 , 0.5-material.x2 ;  %14
                                                       0.5+material.x2 , 0.5-material.x2 ,    -material.x2 ;  %12
                                                          -material.x2 ,    -material.x2 ,    -material.x2 ;  %13
                                                           material.x2 , 0.5-material.x2 , 0.5+material.x2 ] ;%15 
material.parameters{2}.cylinder_top_centers = mod(material.parameters{2}.cylinder_top_centers,1);
material.parameters{2}.cylinder_top_centers = 0.5*(material.parameters{2}.cylinder_bot_centers+material.parameters{2}.cylinder_top_centers);

material.parameters{2}.cylinder_radius = cylinder_radius(2)*ones(1,size(material.parameters{2}.cylinder_bot_centers,1)); 
%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};