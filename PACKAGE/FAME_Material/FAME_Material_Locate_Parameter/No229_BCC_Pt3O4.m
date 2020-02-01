function material = No229_BCC_Pt3O4(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No229_BCC_Pt3O4
%
% This is a simple defect superstructure of the CsCl (B2) structure. One atom has been removed from a 2กักั2กักั2 supercell of CsCl.
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A4B3_cI14_229_c_b.html
%
% Edit at 2017/7/19 By Hsiao-Han Huang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'Pt3O4';
material.material_num = 2;

material.lattice_type = 'body_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A4B3_cI14_229_c_b.html';

material.lattice_constant.a = 6.226; 
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Pt
material.parameters{1}.name                 = 'Pt';
material.parameters{1}.color_map            = [192 192 192]/255;

material.parameters{1}.sphere_centers       = [      0,     .5,     .5;
                                                    .5,      0,     .5;
                                                    .5,     .5,      0;
                                              ];
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 

%% Shape description for O
material.parameters{2}.name                 = 'O';
material.parameters{2}.color_map            = [255 0 0]/255;

material.parameters{2}.sphere_centers       = [     .5,     .5,     .5;
                                                     0,      0,     .5;
                                                     0,     .5,      0;
                                                    .5,      0,      0;
                                              ];
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 
 
%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};