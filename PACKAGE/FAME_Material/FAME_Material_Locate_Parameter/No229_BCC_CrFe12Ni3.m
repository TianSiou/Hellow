function material = No229_BCC_CrFe12Ni3(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No229_BCC_CrFe12Ni3
%
% Austenitic steels are alloys of iron and other metals with an averaged face-centered cubic structure. 
% This model is not meant to represent a real steel, and the selection of atom types for each Wyckoff position is arbitrary. 
% If we set the y3=1/4y3=1/4 then the atoms are on the sites of an fcc lattice.
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/AB12C3_cI32_229_a_h_b.html
%
% Edit at 2017/7/19 By Hsiao-Han Huang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'CrFe12Ni3';
material.material_num = 3;

material.lattice_type = 'body_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB12C3_cI32_229_a_h_b.html';
material.y3 = 0.7625;

material.lattice_constant.a = 7.04;  
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

material.parameters{3}.sphere_centers       = [      2*material.y3-1,          material.y3,           material.y3;
                                                                   0,          material.y3,        -material.y3+1;
                                                                   0,       -material.y3+1,           material.y3;
                                                    -2*material.y3+2,       -material.y3+1,        -material.y3+1;
                                                         material.y3,      2*material.y3-1,           material.y3;
                                                      -material.y3+1,                    0,          material.y3;
                                                         material.y3,                    0,       -material.y3+1;
                                                      -material.y3+1,     -2*material.y3+2,       -material.y3+1;
                                                         material.y3,          material.y3,      2*material.y3-1;
                                                         material.y3,       -material.y3+1,                    0;
                                                      -material.y3+1,          material.y3,                    0;
                                                      -material.y3+1,       -material.y3+1,     -2*material.y3+2;
                                              ];
material.parameters{3}.sphere_radius        = sphere_radius(3)*ones(1,size(material.parameters{3}.sphere_centers,1)); 

%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map, material.parameters{3}.color_map};