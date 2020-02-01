function material = No221_SC_BaHg11(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No221_SC_BaHg11
%
% See detail on following website:
%      http://www.aflowlib.org/CrystalDatabase/AB11_cP36_221_c_agij.html
% Edit at 2017/6/22 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'BaHg11';
material.material_num = 2;

material.lattice_type = 'simple_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB11_cP36_221_c_agij.html';
material.x3           = 0.345 ;
material.y4           = 0.225 ;
material.y5           = 0.115 ;

material.lattice_constant.a = 9.60000;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Ba
material.parameters{1}.name                 = 'Ba';
material.parameters{1}.color_map            = [0 221 0]/255;

material.parameters{1}.sphere_centers       = [      0,     0.5,    0.5 ;
                                                  0.5 ,      0 ,    0.5 ;
                                                  0.5 ,    0.5 ,      0 ;  ];
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1));
    
%% Shape description for Hg
material.parameters{2}.name                 = 'Hg';
material.parameters{2}.color_map            = [175 175 216]/255;

material.parameters{2}.sphere_centers       = [              0 ,               0 ,               0 ;
                                                    material.x3 ,     material.x3 ,      material.x3 ;
                                                       1-material.x3 ,  1-material.x3 ,    material.x3 ;
                                                             1-material.x3 ,     material.x3 ,  1-material.x3 ;
                                                               material.x3 ,   1-material.x3 ,   1-material.x3 ;
                                                            material.x3 , material.x3 , 1-material.x3 ;
                                                            1-material.x3 , 1-material.x3 , 1-material.x3 ; 
                                                            material.x3 , 1-material.x3 , material.x3 ;
                                                            1-material.x3 , material.x3 , material.x3 ;
                                                            0 , material.y4 , material.y4 ;
                                                            0 , material.y4 , 1-material.y4 ;
                                                            0 , 1-material.y4 , material.y4 ;
                                                            0 , 1-material.y4 , 1-material.y4 ; 
                                                             material.y4 , 0 , material.y4 ;
                                                            material.y4 , 0 , 1-material.y4 ;
                                                            1-material.y4 , 0 , material.y4 ;
                                                           1-material.y4 , 0 , 1-material.y4 ;
                                                             material.y4 , material.y4 , 0 ;
                                                            material.y4 , 1-material.y4 , 0 ;
                                                            1-material.y4 , material.y4 , 0 ;
                                                           1-material.y4 , 1-material.y4 , 0 ;
                                                           0.5 , material.y5 , material.y5 ;
                                                            0.5 , material.y5 , 1-material.y5 ;
                                                             0.5 , 1-material.y5 , material.y5 ;
                                                             0.5 , 1-material.y5 , 1-material.y5 ;
                                                            material.y5 , 0.5 , material.y5 ;
                                                            material.y5 , 0.5 , 1-material.y5 ;
                                                            1-material.y5 , 0.5 , material.y5 ;
                                                             1-material.y5 , 0.5 , 1-material.y5 ;
                                                            material.y5 , material.y5 , 0.5 ;
                                                             material.y5 , 1-material.y5 , 0.5 ;
                                                              1-material.y5 , material.y5 , 0.5 ;
                                                              1-material.y5 , 1-material.y5 , 0.5 ;];
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1));  

    
%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};
 
  