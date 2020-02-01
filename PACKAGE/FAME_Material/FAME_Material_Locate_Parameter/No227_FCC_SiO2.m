function material = No227_FCC_SiO2(sphere_radius, cylinder_radius)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No227_FCC_SiO2
%
% See detail on following website:
%      http://www.aflowlib.org/CrystalDatabase/A2B_cF24_227_c_a.html
% Edit at 2017/7/9 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'SiO2';
material.material_num = 2;

material.lattice_type = 'face_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A2B_cF24_227_c_a.html';


material.lattice_constant.a = 7.16600;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Ni
material.parameters{1}.name                 = 'Si';
material.parameters{1}.color_map            = [222 191 150]/255;
material.parameters{1}.sphere_centers       = [      0.125 , 0.125 , 0.125 ;
                                                     0.875 , 0.875 , 0.875 ];
                                                   
material.parameters{1}.sphere_centers = mod(material.parameters{1}.sphere_centers,1)  ;                                                                                    
material.parameters{1}.sphere_radius        =sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1));  
material.parameters{1}.cylinder_bot_centers =  [   0.125 , 0.125 , 0.125 ;
                                                   0.125 , 0.125 , 0.125 ;
                                                   0.125 , 0.125 , 0.125 ;
                                                   0.125 , 0.125 , 0.125 ;
                                                   0.875 , 0.875 , 0.875 ;
                                                   0.875 , 0.875 , 0.875 ;
                                                   0.875 , 0.875 , 0.875 ;
                                                   0.875 , 0.875 , 0.875];
                                               
material.parameters{1}.cylinder_top_centers =  [   0 ,    0 ,    0 ;  
                                                   0 ,    0 ,  0.5 ;  
                                                   0 ,  0.5 ,    0 ;
                                                 0.5 ,    0 ,    0 ;
                                                   1 ,    1 ,    1 ;
                                                   1 ,    1 ,  0.5 ;
                                                   1 ,  0.5 ,    1 ;
                                                 0.5 ,    1 ,    1 ];
% material.parameters{1}.cylinder_top_centers = mod(material.parameters{1}.cylinder_top_centers,1);
material.parameters{1}.cylinder_top_centers = 0.5*(material.parameters{1}.cylinder_top_centers+material.parameters{1}.cylinder_bot_centers) ;     
                                                     
material.parameters{1}.cylinder_radius        = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_bot_centers,1));    
%% Shape description for Ti
material.parameters{2}.name                 = 'O';
material.parameters{2}.color_map            =  [255 0 0]/255;
material.parameters{2}.sphere_centers       = [                 0 ,                0 ,                0 ;
                                                                0 ,                0 ,              0.5 ;
                                                                0 ,              0.5 ,                0 ;
                                                              0.5 ,                0 ,                0 ];
                                                
material.parameters{2}.sphere_centers = mod(material.parameters{2}.sphere_centers,1)  ;                                                                                       
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1));  
material.parameters{2}.cylinder_bot_centers = material.parameters{1}.cylinder_top_centers;
                                               
material.parameters{2}.cylinder_top_centers =  [                0 ,     0 ,     0 ;
                                                                0       0 ,   0.5 ;
                                                                0 ,   0.5 ,     0 ;
                                                              0.5 ,     0 ,     0 ;
                                                                1 ,     1 ,     1 ;
                                                                1 ,     1 ,   0.5 ;
                                                                1 ,   0.5 ,     1 ;
                                                              0.5 ,     1 ,     1 ];

% material.parameters{2}.cylinder_top_centers = mod(material.parameters{2}.cylinder_top_centers,1);                                                     
material.parameters{2}.cylinder_radius        = cylinder_radius(2)*ones(1,size(material.parameters{2}.cylinder_bot_centers,1));       


%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map };
end
    