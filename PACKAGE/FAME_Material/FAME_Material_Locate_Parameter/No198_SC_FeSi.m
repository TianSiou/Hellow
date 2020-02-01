function material = No198_SC_FeSi(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No198_SC_FeSi
%
% When x1=0 and x2=1/2, or x1=1/4 and x2=3/4, this lattice reduces to the rock salt (B1) structure. 
% When x1=-x2=1/8*(sqrt(5)-1), we have an ideal structure where every atom is seven-fold coordinated. 
% See detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/AB_cP8_198_a_a.FeSi.html
% Edit at 2017/6/20 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'FeSi';
material.material_num = 2;

material.lattice_type = 'simple_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB_cP8_198_a_a.FeSi.html';
material.x1           =  (1/8)*(sqrt(5)-1);
material.x2           = -(1/8)*(sqrt(5)-1);

material.lattice_constant.a = 4.48688;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Fe
material.parameters{1}.name                 = 'Fe';
material.parameters{1}.color_map            = [128 64 0]/255;
material.parameters{1}.sphere_centers       = [        material.x1 ,      material.x1 ,        material.x1 ;
                                                   0.5-material.x1 ,    1-material.x1 ,    0.5+material.x1 ;
                                                     1-material.x1 ,  0.5+material.x1 ,    0.5-material.x1 ;
                                                   0.5+material.x1 ,  0.5-material.x1 ,      1-material.x1 ];
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 

material.parameters{1}.cylinder_bot_centers =  [      material.x1 ,      material.x1 ,        material.x1 ;
                                                      material.x1 ,      material.x1 ,        material.x1 ;
                                                      material.x1 ,      material.x1 ,        material.x1 ;
                                                  0.5-material.x1 ,    1-material.x1 ,    0.5+material.x1 ;
                                                  0.5-material.x1 ,    1-material.x1 ,    0.5+material.x1 ;
                                                  0.5-material.x1 ,    1-material.x1 ,    0.5+material.x1 ;
                                                    1-material.x1 ,  0.5+material.x1 ,    0.5-material.x1 ;
                                                    1-material.x1 ,  0.5+material.x1 ,    0.5-material.x1 ;
                                                    1-material.x1 ,  0.5+material.x1 ,    0.5-material.x1 ;
                                                  0.5+material.x1 ,  0.5-material.x1 ,      1-material.x1 ;
                                                  0.5+material.x1 ,  0.5-material.x1 ,      1-material.x1 ;
                                                  0.5+material.x1 ,  0.5-material.x1 ,      1-material.x1 ];
                                                   
material.parameters{1}.cylinder_top_centers = 0.5*[    material.x1+ 0.5-material.x2 ,         material.x1-material.x2 ,    material.x1+ 0.5+material.x2 ;
                                                            material.x1-material.x2 ,     material.x1+0.5+material.x2 ,     material.x1+0.5-material.x2 ;
                                                        material.x1+0.5+material.x2 ,     material.x1+0.5-material.x2 ,         material.x1-material.x2 ;
                                                      0.5-material.x1+1+material.x2 ,     1-material.x1+1+material.x2 ,   0.5+material.x1+1+material.x2 ;
                                                        0.5-material.x1-material.x2 ,   1-material.x1+0.5+material.x2 , 0.5+material.x1+0.5-material.x2 ;
                                                    0.5-material.x1+0.5+material.x2 ,   1-material.x1+0.5-material.x2 ,     0.5+material.x1-material.x2 ;
                                                        1-material.x1+1+material.x2 ,   0.5+material.x1+1+material.x2 ,   0.5-material.x1+1+material.x2 ;
                                                      1-material.x1+0.5-material.x2 ,     0.5+material.x1-material.x2 , 0.5-material.x1+0.5+material.x2 ;
                                                      1-material.x1+0.5+material.x2 , 0.5+material.x1+0.5-material.x2 ,     0.5-material.x1-material.x2 ;
                                                      0.5+material.x1+1+material.x2 ,   0.5-material.x1+1+material.x2 ,     1-material.x1+1+material.x2 ;
                                                    0.5+material.x1+0.5-material.x2 ,     0.5-material.x1-material.x2 ,   1-material.x1+0.5+material.x2 ;
                                                        0.5+material.x1-material.x2 , 0.5-material.x1+0.5+material.x2 ,   1-material.x1+0.5-material.x2 ];
material.parameters{1}.cylinder_radius        = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_bot_centers,1)); 
%% Shape description for Si
material.parameters{2}.name                 = 'Si';
material.parameters{2}.color_map            = [223 193 153]/255;
material.parameters{2}.sphere_centers       = [       1+material.x2,     1+material.x2,       1+material.x2;
                                                    0.5-material.x2,      -material.x2,     0.5+material.x2;
                                                       -material.x2,   0.5+material.x2,     0.5-material.x2;
                                                    0.5+material.x2,   0.5-material.x2,        -material.x2];
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 

material.parameters{2}.cylinder_bot_centers = [     0.5-material.x2,      -material.x2,     0.5+material.x2;
                                                       -material.x2,   0.5+material.x2,     0.5-material.x2;
                                                    0.5+material.x2,   0.5-material.x2,        -material.x2;
                                                      1+material.x2,     1+material.x2,       1+material.x2;
                                                       -material.x2,   0.5+material.x2,     0.5-material.x2;
                                                    0.5+material.x2,   0.5-material.x2,        -material.x2
                                                      1+material.x2,     1+material.x2,       1+material.x2; 
                                                    0.5-material.x2,      -material.x2,     0.5+material.x2;
                                                    0.5+material.x2,   0.5-material.x2,        -material.x2; 
                                                      1+material.x2,     1+material.x2,       1+material.x2;
                                                    0.5-material.x2,      -material.x2,     0.5+material.x2;
                                                       -material.x2,   0.5+material.x2,     0.5-material.x2;
                                                    0.5-material.x2,      -material.x2,     0.5+material.x2;
                                                       -material.x2,   0.5+material.x2,     0.5-material.x2;
                                                    0.5+material.x2,   0.5-material.x2,        -material.x2];  
material.parameters{2}.cylinder_top_centers = [  material.parameters{1}.cylinder_top_centers;
                                                       -material.x2,   0.5+material.x2,     0.5-material.x2;
                                                    0.5+material.x2,   0.5-material.x2,        -material.x2;
                                                    0.5-material.x2,      -material.x2,     0.5+material.x2 ];                                              
material.parameters{2}.cylinder_radius        = cylinder_radius(2)*ones(1,size(material.parameters{2}.cylinder_bot_centers,1));
%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};

    