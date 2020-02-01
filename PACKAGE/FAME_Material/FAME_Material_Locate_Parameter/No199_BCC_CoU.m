function material = No199_BCC_CoU(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No199_BCC_CoU
%
% When x1=1/4 and x2=0, or visa versa, this structure reduces to CsCl (B2) with aB2 = 1/2 a$.
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/AB_cI16_199_a_a.html
%
% Edit at 2017/6/19 By Jia-Wei Lin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'CoU';
material.material_num = 2;

material.lattice_type = 'body_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB_cI16_199_a_a.html';
material.x1           = 0.294;
material.x2           = 0.0347;

material.lattice_constant.a = 6.3557;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for C
material.parameters{1}.name                 = 'Co';
material.parameters{1}.color_map            = [47 153 249]/255;

material.parameters{1}.sphere_centers       = [      2*material.x1,      2*material.x1,      2*material.x1;
                                                               0.5,                  0,  1.5-2*material.x1;
                                                                 0,  1.5-2*material.x1,                0.5;
                                                 1.5-2*material.x1,                0.5,                 0];
                                          
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 
%% Shape description for Fe
material.parameters{2}.name                 = 'U';
material.parameters{2}.color_map            = [251 155 186]/255;

material.parameters{2}.sphere_centers       = [      2*material.x2,      2*material.x2,      2*material.x2;
                                                               0.5,                  0,  0.5-2*material.x2;
                                                                 0,  0.5-2*material.x2,                0.5;
                                                 0.5-2*material.x2,                0.5,                  0];
                                             
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 

%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};