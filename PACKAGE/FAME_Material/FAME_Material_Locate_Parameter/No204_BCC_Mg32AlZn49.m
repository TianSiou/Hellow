function material = No204_BCC_Mg32AlZn49(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No204_BCC_Mg32AlZn49
%
% Most of the sites in this lattice have random occupancy. 
% In particular, according to (Bergman, 1957): The Al¡VI (2a) site is only occupied 80% of the time, the Zn¡VI (24g) site is occupied by Al 19% of the time, 
% the Zn¡VII (24g) site is occupied by Al 43% of the time, and the Zn¡VIII (48h) site is occupied by Al 36% of the time.
% More detail on following website:
%       http://www.aflowlib.org/CrystalDatabase/AB32C48_cI162_204_a_2efg_2gh.html
%
% Edit at 2017/7/19 By À¹±á¬v(Check by Hsiao-Han Huang)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'Mg32AlZn49';
material.material_num = 8;

material.lattice_type = 'body_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB32C48_cI162_204_a_2efg_2gh.html';
material.x2           = 0.8203;
material.x3           = 0.5998;
material.x4           = 0.1836;
material.y5           = 0.2942;
material.z5           = 0.8806;
material.y6           = 0.0908;
material.z6           = 0.8499;
material.y7           = 0.1748;
material.z7           = 0.6993;
material.x8           = 0.686;
material.y8           = 0.0969;
material.z8           = 0.332;

material.lattice_constant.a = 14.16;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Al
material.parameters{1}.name                 = 'Al';
material.parameters{1}.color_map            = [128 128 128]/255;

material.parameters{1}.sphere_centers       = [     0,     0,     0];                                             
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 

%% Shape description for Mg1
material.parameters{2}.name                 = 'Mg1';
material.parameters{2}.color_map            = [0 255 0]/255;

material.parameters{2}.sphere_centers       = [                  .5,     .5+material.x2-1,          material.x2;
                                                                 .5,     .5-material.x2+1,       -material.x2+1;
                                                        material.x2,                   .5,     .5+material.x2-1;
                                                     -material.x2+1,                   .5,     .5-material.x2+1;
                                                   .5+material.x2-1,          material.x2,                   .5;
                                                   .5-material.x2+1,       -material.x2+1,                   .5;
                                             ];
material.parameters{2}.sphere_centers=mod(material.parameters{2}.sphere_centers,1);
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 

%% Shape description for Mg2
material.parameters{3}.name                 = 'Mg2';
material.parameters{3}.color_map            = [0 255 0]/255;

material.parameters{3}.sphere_centers       = [                  .5,     .5+material.x3-1,          material.x3;
                                                                 .5,     .5-material.x3+1,       -material.x3+1;
                                                        material.x3,                   .5,     .5+material.x3-1;
                                                     -material.x3+1,                   .5,     .5-material.x3+1;
                                                   .5+material.x3-1,          material.x3,                   .5;
                                                   .5-material.x3+1,       -material.x3+1,                   .5;
                                             ];
material.parameters{3}.sphere_centers=mod(material.parameters{3}.sphere_centers,1);                                          
material.parameters{3}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{3}.sphere_centers,1)); 

%% Shape description for Mg3
material.parameters{4}.name                 = 'Mg3';
material.parameters{4}.color_map            = [0 255 0]/255;

material.parameters{4}.sphere_centers       = [      2*material.x4,        2*material.x4,        2*material.x4;
                                                     2*material.x4,                    0,                    0;
                                                                 0,        2*material.x4,                    0;
                                                                 0,                    0,        2*material.x4;
                                                  -2*material.x4+1,     -2*material.x4+1,     -2*material.x4+1;
                                                  -2*material.x4+1,                    0,                    0;
                                                                 0,     -2*material.x4+1,                    0;
                                                                 0,                    0,     -2*material.x4+1;
                                             ];
material.parameters{4}.sphere_centers=mod(material.parameters{4}.sphere_centers,1);
radius = 0.03*material.lattice_constant.a;                                               
material.parameters{4}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{4}.sphere_centers,1));

%% Shape description for Mg4
material.parameters{5}.name                 = 'Mg4';
material.parameters{5}.color_map            = [0 255 0]/255;

material.parameters{5}.sphere_centers       = [     material.y5+material.z5-1,                   material.z5,                   material.y5;
                                                      material.z5-material.y5,                   material.z5,                -material.y5+1;
                                                    material.y5-material.z5+1,                -material.z5+1,                   material.y5;
                                                   -material.y5-material.z5+2,                -material.z5+1,                -material.y5+1;
                                                                  material.y5,     material.y5+material.z5-1,                   material.z5;
                                                                  material.y5,     material.y5-material.z5+1,                -material.z5+1;
                                                               -material.y5+1,       material.z5-material.y5,                   material.z5;
                                                               -material.y5+1,    -material.y5-material.z5+2,                -material.z5+1;
                                                                  material.z5,                   material.y5,     material.y5+material.z5-1;
                                                                  material.z5,                -material.y5+1,       material.z5-material.y5;
                                                               -material.z5+1,                   material.y5,     material.y5-material.z5+1;
                                                               -material.z5+1,                -material.y5+1,    -material.y5-material.z5+2;
                                             ];
material.parameters{5}.sphere_centers=mod(material.parameters{5}.sphere_centers,1);                                            
material.parameters{5}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{5}.sphere_centers,1));

%% Shape description for Zn1
material.parameters{6}.name                 = 'Zn1';
material.parameters{6}.color_map            = [0 0 255]/255;

material.parameters{6}.sphere_centers       = [        material.y6+material.z6,                    material.z6,                    material.y6;
                                                       material.z6-material.y6,                    material.z6,                 -material.y6+1;
                                                     material.y6-material.z6+1,                 -material.z6+1,                    material.y6;
                                                    -material.y6-material.z6+1,                 -material.z6+1,                 -material.y6+1;
                                                                   material.y6,        material.y6+material.z6,                    material.z6;
                                                                   material.y6,      material.y6-material.z6+1,                 -material.z6+1;
                                                                -material.y6+1,        material.z6-material.y6,                    material.z6;
                                                                -material.y6+1,     -material.y6-material.z6+1,                 -material.z6+1;
                                                                   material.z6,                    material.y6,        material.y6+material.z6;
                                                                   material.z6,                 -material.y6+1,        material.z6-material.y6;
                                                                -material.z6+1,                    material.y6,      material.y6-material.z6+1;
                                                                -material.z6+1,                 -material.y6+1,     -material.y6-material.z6+1;
                                             ];
material.parameters{6}.sphere_centers=mod(material.parameters{6}.sphere_centers,1);                                           
material.parameters{6}.sphere_radius        = sphere_radius(3)*ones(1,size(material.parameters{6}.sphere_centers,1)); 

%% Shape description for Zn2
material.parameters{7}.name                 = 'Zn2';
material.parameters{7}.color_map            = [0 0 255]/255;

material.parameters{7}.sphere_centers       = [       material.y7+material.z7,                    material.z7,                    material.y7;
                                                      material.z7-material.y7,                    material.z7,                 -material.y7+1;
                                                    material.y7-material.z7+1,                 -material.z7+1,                    material.y7;
                                                   -material.y7-material.z7+1,                 -material.z7+1,                 -material.y7+1;
                                                                  material.y7,        material.y7+material.z7,                    material.z7;
                                                                  material.y7,      material.y7-material.z7+1,                 -material.z7+1;
                                                               -material.y7+1,        material.z7-material.y7,                    material.z7;
                                                               -material.y7+1,     -material.y7-material.z7+1,                 -material.z7+1;
                                                                  material.z7,                    material.y7,        material.y7+material.z7;
                                                                  material.z7,                 -material.y7+1,        material.z7-material.y7;
                                                               -material.z7+1,                    material.y7,      material.y7-material.z7+1;
                                                               -material.z7+1,                 -material.y7+1,     -material.y7-material.z7+1;
                                             ];
material.parameters{7}.sphere_centers=mod(material.parameters{7}.sphere_centers,1);                                            
material.parameters{7}.sphere_radius        = sphere_radius(3)*ones(1,size(material.parameters{7}.sphere_centers,1)); 

%% Shape description for Zn3
material.parameters{8}.name                 = 'Zn3';
material.parameters{8}.color_map            = [0 0 255]/255;

material.parameters{8}.sphere_centers       = [       material.y8+material.z8,      material.z8+material.x8-1,        material.x8+material.y8;
                                                      material.z8-material.y8,      material.z8-material.x8+1,     -material.x8-material.y8+1;
                                                    material.y8-material.z8+1,     -material.z8-material.x8+2,      material.y8-material.x8+1;
                                                   -material.y8-material.z8+1,        material.x8-material.z8,        material.x8-material.y8;
                                                   -material.y8-material.z8+1,     -material.z8-material.x8+2,     -material.x8-material.y8+1;
                                                    material.y8-material.z8+1,        material.x8-material.z8,        material.x8+material.y8;
                                                      material.z8-material.y8,      material.z8+material.x8-1,        material.x8-material.y8;
                                                      material.y8+material.z8,      material.z8-material.x8+1,      material.y8-material.x8+1;
                                                      material.x8+material.y8,        material.y8+material.z8,      material.z8+material.x8-1;
                                                    material.y8-material.x8+1,      material.y8-material.z8+1,     -material.z8-material.x8+2;
                                                      material.x8-material.y8,     -material.y8-material.z8+1,        material.x8-material.z8;
                                                   -material.x8-material.y8+1,        material.z8-material.y8,      material.z8-material.x8+1;
                                                   -material.x8-material.y8+1,     -material.y8-material.z8+1,     -material.z8-material.x8+2;
                                                      material.x8-material.y8,        material.z8-material.y8,      material.z8+material.x8-1;
                                                    material.y8-material.x8+1,        material.y8+material.z8,      material.z8-material.x8+1;
                                                      material.x8+material.y8,      material.y8-material.z8+1,        material.x8-material.z8;
                                                    material.z8+material.x8-1,        material.x8+material.y8,        material.y8+material.z8;
                                                      material.x8-material.z8,        material.x8-material.y8,     -material.y8-material.z8+1;
                                                    material.z8-material.x8+1,     -material.x8-material.y8+1,        material.z8-material.y8;
                                                   -material.z8-material.x8+2,      material.y8-material.x8+1,    1+material.y8-material.z8+1;
                                                   -material.z8-material.x8+2,     -material.x8-material.y8+1,     -material.y8-material.z8+1;
                                                    material.z8-material.x8+1,      material.y8-material.x8+1,        material.y8+material.z8;
                                                      material.x8-material.z8,        material.x8+material.y8,      material.y8-material.z8+1;
                                                    material.z8+material.x8-1,        material.x8-material.y8,        material.z8-material.y8;
                                             ];
material.parameters{8}.sphere_centers=mod(material.parameters{8}.sphere_centers,1);                                              
material.parameters{8}.sphere_radius        = sphere_radius(3)*ones(1,size(material.parameters{8}.sphere_centers,1));

%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map,  material.parameters{3}.color_map, material.parameters{4}.color_map, material.parameters{5}.color_map, material.parameters{6}.color_map, material.parameters{7}.color_map, material.parameters{8}.color_map};