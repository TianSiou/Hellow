function material = No217_BCC_Cu5Zn8(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No217_BCC_Cu5Zn8
%
% See detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A5B8_cI52_217_ce_cg.html
% Edit at 2017/6/22 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'Cu5Zn8';
material.material_num = 2;

material.lattice_type = 'body_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A5B8_cI52_217_ce_cg.html';
material.x1           = 0.32774  ;
material.x2           = 0.10781 ;
material.x3           = 0.64421  ;
material.x4           = 0.68844  ;
material.z4           = 0.03674  ;

material.lattice_constant.a = 8.8664;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Zn
material.parameters{1}.name                 = 'Zn';
material.parameters{1}.color_map            = [105 81 189]/255;
material.parameters{1}.sphere_centers       =  [          2*material.x2 ,               2*material.x2 ,              2*material.x2 ;
                                                                      0 ,                           0 ,            1-2*material.x2 ;
                                                                      0 ,             1-2*material.x2 ,                          0 ;
                                                        1-2*material.x2 ,                           0 ,                          0 ;
                                                material.x4+material.z4 ,     material.x4+material.z4 ,            2*material.x4-1 ;
                                              material.z4-material.x4+1 ,   material.z4-material.x4+1 ,           -2*material.x4+2 ;
                                                material.x4-material.z4 ,  -material.x4-material.z4+1 ,                          0 ;
                                             -material.x4-material.z4+1 ,     material.x4-material.z4 ,                          0 ;
                                                        2*material.x4-1 ,     material.x4+material.z4 ,    material.x4+material.z4 ;
                                                       -2*material.x4+2 ,   material.z4-material.x4+1 ,  material.z4-material.x4+1 ;
                                                                      0 ,     material.x4-material.z4 , -material.x4-material.z4+1 ;
                                                                      0 ,  -material.x4-material.z4+1 ,   material.x4-material.z4  ;
                                                material.x4+material.z4 ,             2*material.x4-1 ,    material.x4+material.z4 ;
                                              material.z4-material.x4+1 ,            -2*material.x4+2 ,  material.z4-material.x4+1 ;
                                             -material.x4-material.z4+1 ,                           0 ,    material.x4-material.z4 ;
                                                material.x4-material.z4 ,                           0 ,  -material.x4-material.z4+1];
                                            
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 
%% Shape description for Cu
material.parameters{2}.name                 = 'Cu';
material.parameters{2}.color_map            = [185 122 87]/255;
material.parameters{2}.sphere_centers       =  [          2*material.x1 ,               2*material.x1 ,              2*material.x1 ;
                                                                      0 ,                           0 ,            1-2*material.x1 ;
                                                                      0 ,             1-2*material.x1 ,                          0 ;
                                                        1-2*material.x1 ,                           0 ,                          0 ;
                                                                      0 ,                 material.x3 ,                material.x3 ;
                                                            material.x3 ,                           0 ,                material.x3 ;
                                                            material.x3 ,                 material.x3 ,                          0 ;
                                                                      0 ,               1-material.x3 ,              1-material.x3 ;
                                                          1-material.x3 ,                           0 ,              1-material.x3 ;
                                                          1-material.x3 ,               1-material.x3 ,                          0 ];                                         
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 
%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};

    