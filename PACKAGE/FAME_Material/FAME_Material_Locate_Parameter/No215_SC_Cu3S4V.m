function material = No215_SC_Cu3S4V(sphere_radius, cylinder_radius)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No215_SC_Cu3S4V
%
%Other compounds with this structure:Cu3S4Nb, Cu3S4Ta, Cu3Se4Nb, Cu3Te4Ta, Cu3Te4V
% See detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A3B4C_cP8_215_d_e_a.html
% Edit at 2017/7/12 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
material.Prototype    = 'Cu3S4V';
material.material_num = 3;

material.lattice_type = 'simple_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A3B4C_cP8_215_d_e_a.html';
material.x3           = 0.2372 ;
material.lattice_constant.a = 5.39120;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for V
material.parameters{1}.name                 = 'V';
material.parameters{1}.color_map            = [127 127 127]/255;
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
material.parameters{2}.sphere_centers       = [  0.5 ,    0 ,     0 ;
                                                   0 ,  0.5 ,     0 ;
                                                   0 ,    0 ,   0.5 ];                                    
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 
material.parameters{2}.cylinder_bot_centers =  [   0.5 ,    0 ,    0 ;
                                                     0 ,  0.5 ,    0 ;
                                                     0 ,    0 ,  0.5 ;
                                                   0.5 ,    1 ,    0 ;
                                                     1 ,  0.5 ,    0 ;
                                                     1 ,    1 ,  0.5 ;
                                                     1 ,  0.5 ,    1 ;
                                                     1 ,    0 ,  0.5 ;
                                                   0.5 ,    0 ,    1 ;
                                                     0 ,    1 ,  0.5 ;
                                                     0 ,  0.5 ,    1 ;
                                                   0.5 ,    1 ,    1 ] ;
                                                                                                   
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
material.parameters{3}.cylinder_top_centers =   [    0 ,    0 ,    0 ;
                                                   0.5 ,    0 ,    0 ;
                                                     0 ,  0.5 ,    0 ;
                                                     0 ,    0 ,  0.5 ;
                                                     1 ,    1 ,    0 ;
                                                   0.5 ,    1 ,    0 ;
                                                     1 ,  0.5 ,    0 ;
                                                     1 ,    1 ,  0.5 ;
                                                     1 ,    0 ,    1 ;
                                                     1 ,  0.5 ,    1 ;
                                                     1 ,    0 ,  0.5 ;
                                                   0.5 ,    0 ,    1 ;
                                                     0 ,    1 ,    1 ;
                                                     0 ,    1 ,  0.5 ;
                                                     0 ,  0.5 ,    1 ;
                                                   0.5 ,    1 ,    1 ] ;
material.parameters{3}.cylinder_top_centers = 0.5*(material.parameters{3}.cylinder_top_centers+material.parameters{3}.cylinder_bot_centers) ;                                                 

material.parameters{3}.cylinder_radius        = cylinder_radius(3)*ones(1,size(material.parameters{3}.cylinder_bot_centers,1)); 



%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map, material.parameters{3}.color_map};