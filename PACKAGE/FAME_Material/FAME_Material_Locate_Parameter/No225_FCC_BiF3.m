function material = No225_FCC_BiF3(sphere_radius, cylinder_radius)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No225_FCC_BiF3
%
% See detail on following website:
%      http://www.aflowlib.org/CrystalDatabase/AB3_cF16_225_a_bc.html
% Edit at 2017/7/11 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'BiF3';
material.material_num = 2;

material.lattice_type = 'face_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB3_cF16_225_a_bc.html';


material.lattice_constant.a = 5.85300;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Bi
material.parameters{1}.name                 = 'Bi';
material.parameters{1}.color_map            = [139 63 139]/255;
material.parameters{1}.sphere_centers       = [  0 , 0 , 0 ];
                                                   
material.parameters{1}.sphere_centers = mod(material.parameters{1}.sphere_centers,1)  ;                                        
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1));  
material.parameters{1}.cylinder_bot_centers =  [      0 ,     0 ,    0 ;
                                                      1 ,     0 ,    0 ;
                                                      0 ,     1 ,    0 ;
                                                      0 ,     0 ,    1 ;
                                                      1 ,     1 ,    1 ;
                                                      0 ,     1 ,    1 ;
                                                      1 ,     0 ,    1 ;
                                                      1 ,     1 ,    0 ];
                                               
material.parameters{1}.cylinder_top_centers =  [     0.25 ,     0.25 ,    0.25 ;
                                                     0.25 ,     0.25 ,    0.25 ;
                                                     0.25 ,     0.25 ,    0.25 ;
                                                     0.25 ,     0.25 ,    0.25 ;
                                                     0.75 ,     0.75 ,    0.75 ;
                                                     0.75 ,     0.75 ,    0.75 ;
                                                     0.75 ,     0.75 ,    0.75 ;
                                                     0.75 ,     0.75 ,    0.75 ];
% material.parameters{1}.cylinder_top_centers = mod(material.parameters{1}.cylinder_top_centers,1);
material.parameters{1}.cylinder_top_centers = 0.5*(material.parameters{1}.cylinder_top_centers+material.parameters{1}.cylinder_bot_centers) ;     
                                                     
material.parameters{1}.cylinder_radius        = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_bot_centers,1));    
%% Shape description for F
material.parameters{2}.name                 = 'F';
material.parameters{2}.color_map            =  [164 210 63]/255;
material.parameters{2}.sphere_centers       = [       0.5 ,      0.5 ,     0.5 ;
                                                     0.25 ,     0.25 ,    0.25 ;
                                                     0.75 ,     0.75 ,    0.75 ];
                                                
material.parameters{2}.sphere_centers = mod(material.parameters{2}.sphere_centers,1)  ;                                            
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1));  
material.parameters{2}.cylinder_bot_centers = [      0.25 ,     0.25 ,    0.25 ;
                                                     0.25 ,     0.25 ,    0.25 ;
                                                     0.25 ,     0.25 ,    0.25 ;
                                                     0.25 ,     0.25 ,    0.25 ;
                                                     0.75 ,     0.75 ,    0.75 ;
                                                     0.75 ,     0.75 ,    0.75 ;
                                                     0.75 ,     0.75 ,    0.75 ;
                                                     0.75 ,     0.75 ,    0.75 ];
                                               
material.parameters{2}.cylinder_top_centers =  material.parameters{1}.cylinder_top_centers ;
material.parameters{2}.cylinder_top_centers = mod(material.parameters{2}.cylinder_top_centers,1);                                                     
material.parameters{2}.cylinder_radius        = cylinder_radius(2)*ones(1,size(material.parameters{2}.cylinder_bot_centers,1));       


%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map };
end
    