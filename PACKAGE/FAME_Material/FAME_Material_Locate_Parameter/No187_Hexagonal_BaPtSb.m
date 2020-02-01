function material = No187_Hexagonal_BaPtSb(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No187_Hexagonal_BaPtSb
%
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/ABC_hP3_187_a_d_f.html
%
% Edit at 2017/7/25 by Jyun-Wei Lin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'BaPtSb';
material.material_num = 3;

material.lattice_type = 'hexagonal';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/ABC_hP3_187_a_d_f.html';

material.lattice_constant.a = 4.535;
material.lattice_constant.c = 1.0769570011*material.lattice_constant.a;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Ba
material.parameters{1}.name                 = 'Ba';
material.parameters{1}.color_map            = [0 255 0]/255;
material.parameters{1}.sphere_centers       = [   0,  0,  0
                                              ];
                                              
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 


%% Shape description for Pt
material.parameters{2}.name                 = 'Pt';
material.parameters{2}.color_map            = [192 192 192]/255;
material.parameters{2}.sphere_centers       = [ 1/3 , 2/3 ,  .5] ;
%material.parameters{2}.sphere_centers=mod(material.parameters{2}.sphere_centers,1);
                                              
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 

%% Shape description for Sb
material.parameters{3}.name                 = 'Sb';
material.parameters{3}.color_map            = [160 60 255]/255;
material.parameters{3}.sphere_centers       = [ 2/3 , 1/3 ,  .5] ;
%material.parameters{2}.sphere_centers=mod(material.parameters{2}.sphere_centers,1);
                                           
material.parameters{3}.sphere_radius        = sphere_radius(3)*ones(1,size(material.parameters{3}.sphere_centers,1)); 

%% Collect color map
material.color_map = {material.parameters{1}.color_map,material.parameters{2}.color_map,material.parameters{3}.color_map};