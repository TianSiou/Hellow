function material = No217_BCC_alpha_Mn(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No217_BCC_alpha_Mn
%
% See detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A_cI58_217_ac2g.html
% Edit at 2017/6/21 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'alpha_Mn';
material.material_num = 1;

material.lattice_type = 'body_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A_cI58_217_ac2g.html';
material.x2           =  0.31787  ;
material.x3           = -0.08958 ;
material.z3           =  0.28194  ;
material.x4           =  0.64294  ;
material.z4           =  0.03457  ;

material.lattice_constant.a = 8.911;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Mn
material.parameters{1}.name                 = 'Mn';
material.parameters{1}.color_map            = [163 73 164]/255;
material.parameters{1}.sphere_centers       =  [                      0 ,                           0 ,                          0 ;
                                                          2*material.x2 ,               2*material.x2 ,              2*material.x2 ;
                                                                      0 ,                           0 ,            1-2*material.x2 ;
                                                                      0 ,             1-2*material.x2 ,                          0 ;
                                                        1-2*material.x2 ,                           0 ,                          0 ;
                                                material.x3+material.z3 ,     material.x3+material.z3 ,            1+2*material.x3 ;
                                                material.z3-material.x3 ,     material.z3-material.x3 ,             -2*material.x3 ;
                                              material.x3-material.z3+1 ,  -material.x3-material.z3+1 ,                          0 ;
                                             -material.x3-material.z3+1 ,   material.x3-material.z3+1 ,                          0 ;
                                                        1+2*material.x3 ,     material.x3+material.z3 ,    material.x3+material.z3 ;
                                                         -2*material.x3 ,     material.z3-material.x3 ,    material.z3-material.x3 ;
                                                                      0 ,   material.x3-material.z3+1 , -material.x3-material.z3+1 ;
                                                                      0 ,  -material.x3-material.z3+1 , material.x3-material.z3+1  ;
                                                material.x3+material.z3 ,             1+2*material.x3 ,    material.x3+material.z3 ;
                                                material.z3-material.x3 ,              -2*material.x3 ,    material.z3-material.x3 ;
                                             -material.x3-material.z3+1 ,                           0 ,  material.x3-material.z3+1 ;
                                              material.x3-material.z3+1 ,                           0 , -material.x3-material.z3+1 ;
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

%% Collect color map
material.color_map = {material.parameters{1}.color_map};

    