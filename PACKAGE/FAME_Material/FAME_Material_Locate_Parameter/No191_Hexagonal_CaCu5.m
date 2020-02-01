function material = No191_Hexagonal_CaCu5(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No191_Hexagonal_CaCu5
%
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/AB5_hP6_191_a_cg.html
%
% Edit at 2017/7/25 by Jyun-Wei Lin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'CaCu5';
material.material_num = 2;

material.lattice_type = 'hexagonal';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB5_hP6_191_a_cg.html';

material.lattice_constant.a = 5.405;
material.lattice_constant.c = 0.773913043478*material.lattice_constant.a;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Ca
material.parameters{1}.name                 = 'Ca';
material.parameters{1}.color_map            = [0 255 0]/255;
material.parameters{1}.sphere_centers       = [   0,  0,  0
                                              ];
                                             
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 


%% Shape description for Cu
material.parameters{2}.name                 = 'Cu';
material.parameters{2}.color_map            = [128 64 0]/255;
material.parameters{2}.sphere_centers       = [ 1/3 , 2/3 ,  0;
                                                2/3 , 1/3 ,  0;
                                                0.5 ,   0 , .5;
                                                  0 ,  .5 , .5;
                                                 .5 ,  .5 , .5]; 
%material.parameters{2}.sphere_centers=mod(material.parameters{2}.sphere_centers,1);
                                            
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 

%% Collect color map
material.color_map = {material.parameters{1}.color_map,material.parameters{2}.color_map};