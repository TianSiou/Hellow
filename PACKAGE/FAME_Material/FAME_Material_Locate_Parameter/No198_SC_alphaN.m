function material = No198_SC_alphaN(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No198_SC_alphaN
%
% See detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A_cP8_198_2a.html
%
% Edit at 2017/6/16 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'alphaN';
material.material_num = 1;

material.lattice_type = 'simple_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A_cP8_198_2a.html';
material.x1           =  0.0699;
material.x2           = -0.0378;

material.lattice_constant.a = 5.65;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Fe
material.parameters{1}.name                 = 'N';
material.parameters{1}.color_map            = [37 37 197]/255;
material.parameters{1}.sphere_centers       = [        material.x1 ,      material.x1 ,        material.x1 ;
                                                   0.5-material.x1 ,    1-material.x1 ,    0.5+material.x1 ;
                                                     1-material.x1 ,  0.5+material.x1 ,    0.5-material.x1 ;
                                                   0.5+material.x1 ,  0.5-material.x1 ,      1-material.x1 ;
                                                     1+material.x2 ,    1+material.x2 ,      1+material.x2 ;
                                                   0.5-material.x2 ,    1-material.x2 ,    0.5+material.x2 ;
                                                     1-material.x2 ,  0.5+material.x2 ,    0.5-material.x2 ;
                                                   0.5+material.x2 ,  0.5-material.x2 ,      1-material.x2];
                                             
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 

%% Collect color map
material.color_map = {material.parameters{1}.color_map};
