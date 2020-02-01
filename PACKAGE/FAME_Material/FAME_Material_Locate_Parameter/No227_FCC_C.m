function material = No227_FCC_C(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No227_FCC_C
%
% This is the first crystal structure to be determined by X-ray diffraction. 
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A_cF8_227_a.html
%
% Edit at 2017/7/16 By Hsiao-Han Huang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'C';
material.material_num = 1;

material.lattice_type = 'face_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A_cF8_227_a.html';

material.lattice_constant.a = 3.55;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for C
material.parameters{1}.name                 = 'C';
material.parameters{1}.color_map            = [124 124 124]/255;
material.parameters{1}.sphere_centers       = [   1/8,   1/8,   1/8; 
                                                  7/8,   7/8,   7/8];                                              
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 

material.parameters{1}.cylinder_bot_centers = [   1/8,   1/8+1,   1/8;
                                                  1/8,   1/8+1,   1/8;
                                                  1/8+1,   1/8,   1/8;
                                                  1/8+1,   1/8,   1/8;
                                                  1/8,   1/8,   1/8+1;
                                                  1/8,   1/8,   1/8+1;
                                                  ];
material.parameters{1}.cylinder_top_centers = [   7/8-1,   7/8,   7/8; 
                                                  7/8,   7/8,   7/8-1;
                                                  7/8,   7/8,   7/8-1;
                                                  7/8,   7/8-1,   7/8;
                                                  7/8,   7/8-1,   7/8;
                                                  7/8-1,   7/8,   7/8; 
                                                   ];                                               
material.parameters{1}.cylinder_radius        = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_bot_centers,1));  
%% Collect color map
material.color_map = {material.parameters{1}.color_map};