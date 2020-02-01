function material = No220_BCC_Li(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No220_BCC_Li
%
% See detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A_cI16_220_c.html
% Edit at 2017/6/22 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'Li';
material.material_num = 1;

material.lattice_type = 'body_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A_cI16_220_c.html';
material.x1           = 0.049  ;

material.lattice_constant.a = 5.2716;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Mn
material.parameters{1}.name                 = 'Li';
material.parameters{1}.color_map            = [134 115 202]/255;
material.parameters{1}.sphere_centers       =  [      2*material.x1 ,      2*material.x1 ,     2*material.x1 ;
                                                                0.5 ,                  0 , 0.5-2*material.x1 ;
                                                                  0 ,  0.5-2*material.x1 ,               0.5 ;
                                                  0.5-2*material.x1 ,                0.5 ,                 0 ;
                                                  0.5+2*material.x1 ,  0.5+2*material.x1 , 0.5+2*material.x1 ; 
                                                                0.5 ,                  0 ,   1-2*material.x1 ; 
                                                    1-2*material.x1 ,                0.5 ,                 0 ;
                                                                  0 ,    1-2*material.x1 ,               0.5 ];                                             
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 

%% Collect color map
material.color_map = {material.parameters{1}.color_map};

    