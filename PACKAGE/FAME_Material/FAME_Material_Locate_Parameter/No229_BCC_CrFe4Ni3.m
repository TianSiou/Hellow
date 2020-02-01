function material = No229_BCC_CrFe4Ni3(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No229_BCC_CrFe4Ni3
%
% This structure represents one possible ordering which might be found in an Fe¡VNi¡VCr steel. 
% If we replace the Cr atoms by Ni, this becomes the CsCl (B2) structure. 
% If we replace both the Cr and Ni atoms by Fe, we get the bcc (A2) structure.
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/AB4C3_cI16_229_a_c_b.html
%
% Edit at 2017/7/19 By Hsiao-Han Huang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'CrFe4Ni3';
material.material_num = 3;

material.lattice_type = 'body_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB4C3_cI16_229_a_c_b.html';

material.lattice_constant.a = 5.74;  
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Cr
material.parameters{1}.name                 = 'Cr';
material.parameters{1}.color_map            = [74 165 255]/255;

material.parameters{1}.sphere_centers       = [      0,     0,     0;
                                              ];
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 

%% Shape description for Ni
material.parameters{2}.name                 = 'Ni';
material.parameters{2}.color_map            = [0 166 0]/255;

material.parameters{2}.sphere_centers       = [      0,     .5,     .5;
                                                    .5,      0,     .5;
                                                    .5,     .5,      0;
                                              ];
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 
 
%% Shape description for Fe
material.parameters{3}.name                 = 'Fe';
material.parameters{3}.color_map            = [255 128 63]/255;

material.parameters{3}.sphere_centers       = [     .5,     .5,     .5;
                                                     0,      0,     .5;
                                                     0,     .5,      0;
                                                    .5,      0,      0;
                                              ];
material.parameters{3}.sphere_radius        = sphere_radius(3)*ones(1,size(material.parameters{3}.sphere_centers,1)); 

%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map, material.parameters{3}.color_map};