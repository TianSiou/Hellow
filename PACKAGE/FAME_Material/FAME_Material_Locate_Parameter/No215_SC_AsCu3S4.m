function material = No215_SC_AsCu3S4(sphere_radius, cylinder_radius)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No215_SC_AsCu3S4
%
% See detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/AB3C4_cP8_215_a_c_e.html
% Edit at 2017/7/11 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
material.Prototype    = 'AsCu3S4';
material.material_num = 3;

material.lattice_type = 'simple_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB3C4_cP8_215_a_c_e.html';
material.x3           = 0.25 ;
material.lattice_constant.a = 5.28000;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for As
material.parameters{1}.name                 = 'As';
material.parameters{1}.color_map            = [181 89 181]/255;
material.parameters{1}.sphere_centers       = [    0 ,   0 ,    0 ];
                                                  

material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 
material.parameters{1}.cylinder_bot_centers =  [ 0 ,   0 ,   0 ;
                                                 1 ,   1 ,   0 ;
                                                 1 ,   0 ,   1 ;
                                                 0 ,   1 ,   1 ];
                                               
material.parameters{1}.cylinder_top_centers =  [   material.x3 ,   material.x3 ,    material.x3 ;
                                                  -material.x3 ,  -material.x3 ,    material.x3 ;
                                                  -material.x3 ,   material.x3 ,   -material.x3 ;
                                                   material.x3 ,  -material.x3 ,   -material.x3 ];
material.parameters{1}.cylinder_top_centers = mod(material.parameters{1}.cylinder_top_centers,1);
material.parameters{1}.cylinder_top_centers = 0.5*(material.parameters{1}.cylinder_top_centers+material.parameters{1}.cylinder_bot_centers);

material.parameters{1}.cylinder_radius        = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_bot_centers,1)); 


%% Shape description for Cu
material.parameters{2}.name                 = 'Cu';
material.parameters{2}.color_map            = [164 102 68]/255;
material.parameters{2}.sphere_centers       = [    0 ,    0.5 ,     0.5 ;
                                                 0.5 ,      0 ,     0.5 ;
                                                 0.5 ,    0.5 ,       0 ];                                    

material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 
material.parameters{2}.cylinder_bot_centers =  [   0.5 ,  0.5 ,    0 ;
                                                     0 ,  0.5 ,  0.5 ;
                                                   0.5 ,    0 ,  0.5 ;
                                                   0.5 ,  0.5 ,    0 ;
                                                     1 ,  0.5 ,  0.5 ;
                                                   0.5 ,    1 ,  0.5 ;
                                                     1 ,  0.5 ,  0.5 ;
                                                   0.5 ,    0 ,  0.5 ;
                                                   0.5 ,  0.5 ,    1 ;
                                                     0 ,  0.5 ,  0.5 ;
                                                   0.5 ,    1 ,  0.5 ;
                                                   0.5 ,  0.5 ,    1 ] ;
                                                                                                   
material.parameters{2}.cylinder_top_centers =  [   material.x3 ,   material.x3 ,    material.x3 ;
                                                   material.x3 ,   material.x3 ,    material.x3 ;
                                                   material.x3 ,   material.x3 ,    material.x3 ;
                                                  -material.x3 ,  -material.x3 ,    material.x3 ;
                                                  -material.x3 ,  -material.x3 ,    material.x3 ;
                                                  -material.x3 ,  -material.x3 ,    material.x3 ;
                                                  -material.x3 ,   material.x3 ,   -material.x3 ;
                                                  -material.x3 ,   material.x3 ,   -material.x3 ;
                                                  -material.x3 ,   material.x3 ,   -material.x3 ;
                                                   material.x3 ,  -material.x3 ,   -material.x3 ;
                                                   material.x3 ,  -material.x3 ,   -material.x3 ;
                                                   material.x3 ,  -material.x3 ,   -material.x3 ];
material.parameters{2}.cylinder_top_centers = mod(material.parameters{2}.cylinder_top_centers,1);
material.parameters{2}.cylinder_top_centers = 0.5*(material.parameters{2}.cylinder_top_centers+material.parameters{2}.cylinder_bot_centers);
radius = 0.05*material.lattice_constant.a;
material.parameters{2}.cylinder_radius = cylinder_radius(2)*ones(1,size(material.parameters{2}.cylinder_bot_centers,1)); 
%% Shape description for S
material.parameters{3}.name                 = 'S';
material.parameters{3}.color_map            = [255 242 0]/255;
material.parameters{3}.sphere_centers       = [    material.x3 ,   material.x3 ,    material.x3 ;
                                                  -material.x3 ,  -material.x3 ,    material.x3 ;
                                                  -material.x3 ,   material.x3 ,   -material.x3 ;
                                                   material.x3 ,  -material.x3 ,   -material.x3 ];
material.parameters{3}.sphere_centers = mod(material.parameters{3}.sphere_centers,1);                                                   

material.parameters{3}.sphere_radius        = sphere_radius(3)*ones(1,size(material.parameters{3}.sphere_centers,1)); 
material.parameters{3}.cylinder_bot_centers =  [  material.x3 ,   material.x3 ,    material.x3 ;
                                                  material.x3 ,   material.x3 ,    material.x3 ;
                                                  material.x3 ,   material.x3 ,    material.x3 ;
                                                  material.x3 ,   material.x3 ,    material.x3 ;
                                                 -material.x3 ,  -material.x3 ,    material.x3 ;
                                                 -material.x3 ,  -material.x3 ,    material.x3 ;
                                                 -material.x3 ,  -material.x3 ,    material.x3 ;
                                                 -material.x3 ,  -material.x3 ,    material.x3 ;
                                                 -material.x3 ,   material.x3 ,   -material.x3 ;
                                                 -material.x3 ,   material.x3 ,   -material.x3 ;
                                                 -material.x3 ,   material.x3 ,   -material.x3 ;
                                                 -material.x3 ,   material.x3 ,   -material.x3 ;
                                                  material.x3 ,  -material.x3 ,   -material.x3 ;
                                                  material.x3 ,  -material.x3 ,   -material.x3 ;
                                                  material.x3 ,  -material.x3 ,   -material.x3 ;
                                                  material.x3 ,  -material.x3 ,   -material.x3 ];
material.parameters{3}.cylinder_bot_centers = mod(material.parameters{3}.cylinder_bot_centers,1);                                               
material.parameters{3}.cylinder_top_centers =  [     0 ,    0 ,   0 ;
                                                   0.5 ,  0.5 ,   0 ;
                                                     0 ,  0.5 , 0.5 ;
                                                   0.5 ,    0 , 0.5 ;
                                                     1 ,   1  ,   0 ;
                                                   0.5 ,  0.5 ,   0 ;
                                                     1 ,  0.5 , 0.5 ;
                                                   0.5 ,    1 , 0.5 ;
                                                     1 ,    0 ,   1 ;
                                                     1 ,  0.5 , 0.5 ;
                                                   0.5 ,    0 , 0.5 ;
                                                   0.5 ,  0.5 ,   1 ;
                                                     0 ,    1 ,   1 ;
                                                     0 ,  0.5 , 0.5 ;
                                                   0.5 ,    1 , 0.5 ;
                                                   0.5 ,  0.5 ,   1 ];
                     
material.parameters{3}.cylinder_top_centers = 0.5*(material.parameters{3}.cylinder_top_centers+material.parameters{3}.cylinder_bot_centers) ;                                                 

material.parameters{3}.cylinder_radius        = cylinder_radius(3)*ones(1,size(material.parameters{3}.cylinder_bot_centers,1)); 



%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map, material.parameters{3}.color_map};