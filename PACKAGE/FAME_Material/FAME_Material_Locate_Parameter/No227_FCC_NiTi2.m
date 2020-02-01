function material = No227_FCC_NiTi2(sphere_radius, cylinder_radius)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No227_FCC_NiTi2
%
% See detail on following website:
%      http://www.aflowlib.org/CrystalDatabase/AB2_cF96_227_e_cf.html
% Edit at 2017/6/28 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'NiTi2';
material.material_num = 2;

material.lattice_type = 'face_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB2_cF96_227_e_cf.html';
material.x2           = 0.215 ;
material.x3           = 0.44 ;

material.lattice_constant.a = 11.27800;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Ni
material.parameters{1}.name                 = 'Ni';
material.parameters{1}.color_map            = [0 255 0]/255;
material.parameters{1}.sphere_centers       = [       material.x2 ,       material.x2 ,       material.x2 ;
                                                      material.x2 ,       material.x2 , 0.5-3*material.x2 ;
                                                      material.x2 , 0.5-3*material.x2 ,       material.x2 ;
                                                0.5-3*material.x2 ,       material.x2 ,       material.x2 ;
                                                     -material.x2 ,      -material.x2 , 0.5+3*material.x2 ;
                                                     -material.x2 ,      -material.x2 ,       -material.x2;
                                                     -material.x2 , 0.5+3*material.x2 ,      -material.x2 ;
                                                0.5+3*material.x2 ,      -material.x2 ,      -material.x2 ];
material.parameters{1}.sphere_centers = mod(material.parameters{1}.sphere_centers,1)  ;                                        
material.parameters{1}.sphere_radius        =sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1));  
   
%% Shape description for Ti
material.parameters{2}.name                 = 'Ti';
material.parameters{2}.color_map            =  [195 195 195]/255;
material.parameters{2}.sphere_centers       = [                 0 ,                0 ,                0 ;
                                                                0 ,                0 ,              0.5 ;
                                                                0 ,              0.5 ,                0 ;
                                                              0.5 ,                0 ,                0 ;
                                                 0.25-material.x3 ,      material.x3 ,      material.x3 ;
                                                      material.x3 , 0.25-material.x3 , 0.25-material.x3 ;
                                                      material.x3 , 0.25-material.x3 ,      material.x3 ;
                                                 0.25-material.x3 ,      material.x3 , 0.25-material.x3 ;
                                                      material.x3 ,      material.x3 , 0.25-material.x3 ;
                                                 0.25-material.x3 , 0.25-material.x3 ,      material.x3 ;
                                                 material.x3+0.75 ,     -material.x3 , material.x3+0.75 ;
                                                     -material.x3 , material.x3+0.75 ,     -material.x3 ; 
                                                     -material.x3 , material.x3+0.75 , material.x3+0.75 ; 
                                                 material.x3+0.75 ,     -material.x3 ,     -material.x3 ;
                                                     -material.x3 ,     -material.x3 , material.x3+0.75 ;
                                                 material.x3+0.75 , material.x3+0.75 ,     -material.x3 ];
material.parameters{2}.sphere_centers = mod(material.parameters{2}.sphere_centers,1)  ;                                                                                        
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1));  
   


%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map };
end
    