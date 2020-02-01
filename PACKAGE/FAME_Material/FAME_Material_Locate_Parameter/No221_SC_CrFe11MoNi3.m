function material = No221_SC_CrFe11MoNi3(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No221_SC_CrFe11MoNi3
%
% See detail on following website:
%      http://www.aflowlib.org/CrystalDatabase/AB11CD3_cP16_221_a_dg_b_c.html
% Edit at 2017/6/22 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype     = 'CrFe11MoNi3';
material.material_num = 4;

material.lattice_type = 'simple_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB11CD3_cP16_221_a_dg_b_c.html';
material.x5           = 0.245 ;
material.lattice_constant.a = 5.74000;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Cr
material.parameters{1}.name                 = 'Cr';
material.parameters{1}.color_map            = [104 104 255]/255;
material.parameters{1}.sphere_centers       = [0 , 0 , 0];
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1));
    
    
%% Shape description for Fe
material.parameters{2}.name                 = 'Fe';
material.parameters{2}.color_map            = [188 96 69]/255;

material.parameters{2}.sphere_centers       = [              0.5 ,              0 ,             0 ;
                                                               0 ,            0.5 ,             0 ;
                                                               0 ,              0 ,           0.5 ;
                                                     material.x5 ,    material.x5 ,   material.x5 ;
                                                   1-material.x5 ,  1-material.x5 ,   material.x5 ;
                                                   1-material.x5 ,    material.x5 , 1-material.x5 ;
                                                     material.x5 ,  1-material.x5 , 1-material.x5 ;
                                                     material.x5 ,    material.x5 , 1-material.x5 ;
                                                   1-material.x5 ,  1-material.x5 , 1-material.x5 ; 
                                                     material.x5 ,  1-material.x5 ,   material.x5 ;
                                                   1-material.x5 ,    material.x5 ,   material.x5 ];
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1));
                                                
%% Shape description for Mo
material.parameters{3}.name                 = 'Mo';
material.parameters{3}.color_map            = [0 162 232]/255;
material.parameters{3}.sphere_centers       = [0.5 , 0.5 , 0.5];
material.parameters{3}.sphere_radius        = sphere_radius(3)*ones(1,size(material.parameters{3}.sphere_centers,1));
        
    
%% Shape description for Ni
material.parameters{4}.name                 = 'Ni';
material.parameters{4}.color_map            = [0 232 0]/255;

material.parameters{4}.sphere_centers       = [       0 ,      0.5 ,     0.5 ;
                                                    0.5 ,       0 ,      0.5 ;
                                                    0.5 ,      0.5 ,       0  ];
material.parameters{4}.sphere_radius        = sphere_radius(4)*ones(1,size(material.parameters{4}.sphere_centers,1));
        
    
%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map, material.parameters{3}.color_map, material.parameters{4}.color_map};
 
  