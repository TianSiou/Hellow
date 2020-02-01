function material = No229_BCC_Sb2Tl7(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No229_BCC_Sb2Tl7
%
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A2B7_cI54_229_e_afh.html
%
% Edit at 2017/7/19 By Hsiao-Han Huang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'Sb2Tl7';
material.material_num = 4;

material.lattice_type = 'body_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A2B7_cI54_229_e_afh.html';
material.x2 = 0.6862;
material.x3 = 0.1704;
material.y4 = 0.6503;

material.lattice_constant.a = 11.618; 
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Tl1
material.parameters{1}.name                 = 'Tl1';
material.parameters{1}.color_map            = [103 52 52]/255;

material.parameters{1}.sphere_centers       = [      0,     0,     0;
                                              ];                                            
material.parameters{1}.sphere_radius        =sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 

%% Shape description for Sb
material.parameters{2}.name                 = 'Sb';
material.parameters{2}.color_map            = [128 0 255]/255;

material.parameters{2}.sphere_centers       = [               0,        material.x2,        material.x2;
                                                    material.x2,                  0,        material.x2;
                                                    material.x2,        material.x2,                  0;
                                                              0,     -material.x2+1,     -material.x2+1;
                                                 -material.x2+1,                  0,     -material.x2+1;
                                                 -material.x2+1,     -material.x2+1,                  0;
                                              ];
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 
 
%% Shape description for Tl2
material.parameters{3}.name                 = 'Tl2';
material.parameters{3}.color_map            = [103 52 52]/255;

material.parameters{3}.sphere_centers       = [        2*material.x3,           2*material.x3,           2*material.x3;
                                                                   0,                       0,        -2*material.x3+1;
                                                                   0,        -2*material.x3+1,                       0;
                                                    -2*material.x3+1,                       0,                       0;
                                                                   0,                       0,           2*material.x3;
                                                    -2*material.x3+1,        -2*material.x3+1,        -2*material.x3+1;
                                                                   0,           2*material.x3,                       0;
                                                       2*material.x3,                       0,                       0;
                                              ];                                              
material.parameters{3}.sphere_radius        = sphere_radius(3)*ones(1,size(material.parameters{3}.sphere_centers,1)); 

%% Shape description for Tl3
material.parameters{4}.name                 = 'Tl3';
material.parameters{4}.color_map            = [103 52 52]/255;

material.parameters{4}.sphere_centers       = [      2*material.y4-1,          material.y4,          material.y4;
                                                                   0,          material.y4,       -material.y4+1;
                                                                   0,       -material.y4+1,          material.y4;
                                                    -2*material.y4+2,       -material.y4+1,       -material.y4+1;
                                                         material.y4,      2*material.y4-1,          material.y4;
                                                      -material.y4+1,                    0,          material.y4;
                                                         material.y4,                    0,       -material.y4+1;
                                                      -material.y4+1,     -2*material.y4+2,       -material.y4+1;
                                                         material.y4,          material.y4,      2*material.y4-2;
                                                         material.y4,       -material.y4+1,                    0;
                                                      -material.y4+1,          material.y4,                    0;
                                                      -material.y4+1,       -material.y4+1,     -2*material.y4+2;
                                              ];
material.parameters{4}.sphere_radius        = sphere_radius(4)*ones(1,size(material.parameters{4}.sphere_centers,1)); 

%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map, material.parameters{3}.color_map, material.parameters{4}.color_map};