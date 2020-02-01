function material = No217_BCC_SiF4(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No217_BCC_SiF4
%
% We determined the lattice constant for this structure from the internal coordinates and the Si¡VF bond length given in the reference.
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A4B_cI10_217_c_a.html
%
% Edit at 2017/7/16 By Hsiao-Han Huang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'SiF4';
material.material_num = 2;

material.lattice_type = 'body_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A4B_cI10_217_c_a.html';
material.x2           = 0.165;

material.lattice_constant.a = 5.45858;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Si
material.parameters{1}.name                 = 'Si';
material.parameters{1}.color_map            = [172 89 89]/255;

material.parameters{1}.sphere_centers       = [     0,     0,     0];                                           
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 

material.parameters{1}.cylinder_bot_centers = [   0,     0,     0;
                                                  1,     0,     0;
                                                  0,     1,     0;
                                                  0,     0,     1;
                                                      ];
                                              
material.parameters{1}.cylinder_top_centers = [   material.x2,      material.x2,      material.x2;
                                               -material.x2+1,                0,                0;
                                                            0,   -material.x2+1,                0; 
                                                            0,                0,   -material.x2+1;
                                                  ];
                                             
material.parameters{1}.cylinder_radius        =cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_bot_centers,1)); 
%% Shape description for F
material.parameters{2}.name                 = 'F';
material.parameters{2}.color_map            = [22 177 5]/255;

material.parameters{2}.sphere_centers       = [      2*material.x2,      2*material.x2,      2*material.x2;
                                                                 0,                  0,   -2*material.x2+1;
                                                                 0,   -2*material.x2+1,                  0;
                                                  -2*material.x2+1,                  0,                  0];                                              
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 

material.parameters{2}.cylinder_bot_centers = [      2*material.x2,      2*material.x2,      2*material.x2;
                                                  -2*material.x2+1,                  0,                  0;
                                                                 0,   -2*material.x2+1,                  0;
                                                                 0,                  0,   -2*material.x2+1;
                                                 ];
                                              
material.parameters{2}.cylinder_top_centers = [      material.x2,      material.x2,      material.x2;
                                                  -material.x2+1,                0,                0;
                                                               0,   -material.x2+1,                0;
                                                               0,                0,   -material.x2+1;
                                                   ];
                                        
material.parameters{2}.cylinder_radius        = cylinder_radius(2)*ones(1,size(material.parameters{2}.cylinder_bot_centers,1)); 
%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};