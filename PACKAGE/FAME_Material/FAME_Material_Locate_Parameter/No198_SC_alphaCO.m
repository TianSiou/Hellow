function material = No198_SC_alphaCO(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No198_SC_alphaCO
%
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/AB_cP8_198_a_a.alpha-CO.html
%
% Edit at 2017/6/16 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'alpha-CO';
material.material_num = 2;

material.lattice_type = 'simple_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB_cP8_198_a_a.alpha-CO.html';
material.x1           = -0.042;
material.x2           =  0.067;

material.lattice_constant.a = 5.63;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for C
material.parameters{1}.name                 = 'C';
material.parameters{1}.color_map            = [185 175 166]/255;
material.parameters{1}.sphere_centers       = [   1+material.x1,   1+material.x1,   1+material.x1;
                                                0.5-material.x1,    -material.x1, 0.5+material.x1;
                                                   -material.x1, 0.5+material.x1, 0.5-material.x1;
                                                0.5+material.x1, 0.5-material.x1,    -material.x1];
                                          
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 
%% Shape description for Fe
material.parameters{2}.name                 = 'O';
material.parameters{2}.color_map            = [244 67 39]/255;
material.parameters{2}.sphere_centers       = [     material.x2,     material.x2,     material.x2;
                                                0.5-material.x2,   1-material.x2, 0.5+material.x2;
                                                  1-material.x2, 0.5+material.x2, 0.5-material.x2;
                                                0.5+material.x2, 0.5-material.x2,   1-material.x2];
                                          
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 

%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};