function material = No197_BCC_Ga4Ni(sphere_radius, cylinder_radius)    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No197_BCC_Ga4Ni
%
% See detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A4B_cI40_197_cde_c.html
% Edit at 2017/6/20 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'Ga4Ni';
material.material_num = 2;

material.lattice_type = 'body_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A4B_cI40_197_cde_c.html';
material.x1           = 0.1668 ;
material.x2           = 0.3345 ;
material.x3           = 0.6476 ;
material.x4           = 0.7484 ;

material.lattice_constant.a = 8.4295;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Ga
material.parameters{1}.name                 = 'Ga';
material.parameters{1}.color_map            = [177 99 99]/255;
material.parameters{1}.sphere_centers       = [       2*material.x1 ,     2*material.x1 ,        2*material.x1 ;
                                                                  0 ,                 0 ,      1-2*material.x1 ;
                                                                  0 ,   1-2*material.x1 ,                    0 ; 
                                                    1-2*material.x1 ,                 0 ,                    0 ;
                                                                  0 ,       material.x3 ,          material.x3 ;
                                                        material.x3 ,                 0 ,          material.x3 ;
                                                        material.x3 ,       material.x3 ,                    0 ;
                                                                  0 ,     1-material.x3 ,        1-material.x3 ;
                                                      1-material.x3 ,                 0 ,        1-material.x3 ;
                                                      1-material.x3 ,     1-material.x3 ,                    0 ;
                                                                0.5 ,       material.x4 ,    0.5+material.x4-1 ;
                                                  0.5+material.x4-1 ,               0.5 ,          material.x4 ;
                                                        material.x4 , 0.5+material.x4-1 ,                  0.5 ;
                                                                0.5 ,    1-material.x4  ,    0.5-material.x4+1 ;
                                                  0.5-material.x4+1 ,               0.5 ,        1-material.x4 ;
                                                      1-material.x4 , 0.5-material.x4+1 ,                  0.5 ];
                                              
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 
%% Shape description for H
material.parameters{2}.name                 = 'Ni';
material.parameters{2}.color_map            = [90 240 120]/255;
material.parameters{2}.sphere_centers       =  [    2*material.x2 ,     2*material.x2 ,        2*material.x2 ;
                                                                0 ,                 0 ,      1-2*material.x2 ;
                                                                0 ,   1-2*material.x2 ,                    0 ; 
                                                  1-2*material.x2 ,                 0 ,                    0 ];
                                              
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 

%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};

    