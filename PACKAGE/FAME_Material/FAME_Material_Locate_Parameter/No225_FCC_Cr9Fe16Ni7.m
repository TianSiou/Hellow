function material = No225_FCC_Cr9Fe16Ni7(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No225_FCC_Cr9Fe16Ni7
%
% See detail on following website:
%       http://www.aflowlib.org/CrystalDatabase/A9B16C7_cF128_225_acd_2f_be.html
% Edit at 2017/6/20 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'Cr9Fe16Ni7';
material.material_num = 3;

material.lattice_type = 'face_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A9B16C7_cF128_225_acd_2f_be.html';
material.x5           =  0.25 ;
material.x6           = 0.875 ;
material.x7           = 0.625 ;
material.lattice_constant.a = 11.48000;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for N
material.parameters{1}.name                 = 'Cr';
material.parameters{1}.color_map            = [104 104 255]/255;
material.parameters{1}.sphere_centers       = [        0 ,      0 ,        0 ;
                                                    0.25 ,   0.25 ,     0.25 ;
                                                    0.75 ,   0.75 ,     0.75 ;
                                                     0.5 ,      0 ,        0 ;
                                                       0 ,    0.5 ,      0.5 ;
                                                       0 ,    0.5 ,        0 ;
                                                     0.5 ,      0 ,      0.5 ;
                                                       0 ,      0 ,      0.5 ;
                                                     0.5 ,    0.5 ,        0 ];
                                                
                                           
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1));
    
%% Shape description for H
material.parameters{2}.name                 = 'Fe';
material.parameters{2}.color_map            = [188 96 69]/255;
material.parameters{2}.sphere_centers       = [       material.x6  ,     material.x6  ,        material.x6  ; 
                                                      material.x6  ,      material.x6 ,   -3*material.x6+3  ; 
                                                      material.x6  , -3*material.x6+3 ,        material.x6  ; 
                                                  -3*material.x6+3 ,     material.x6  ,        material.x6  ; 
                                                   -material.x6+1  ,   -material.x6+1 ,    3*material.x6-2  ; 
                                                    -material.x6+1 ,   -material.x6+1 ,      -material.x6+1 ;
                                                   -material.x6+1  ,  3*material.x6-2 ,     -material.x6+1  ; 
                                                   3*material.x6-2 ,  -material.x6+1  ,     -material.x6+1  ; 
                                                      material.x7  ,     material.x7  ,        material.x7  ; 
                                                      material.x7  ,      material.x7 ,   -3*material.x7+2  ; 
                                                      material.x7  , -3*material.x7+2 ,        material.x7  ; 
                                                  -3*material.x7+2 ,      material.x7 ,        material.x7  ; 
                                                   -material.x7+1  ,   -material.x7+1 ,    3*material.x7-1  ; 
                                                    -material.x7+1 ,   -material.x7+1 ,      -material.x7+1 ;
                                                   -material.x7+1  ,  3*material.x7-1 ,     -material.x7+1  ; 
                                                   3*material.x7-1 ,  -material.x7+1  ,      -material.x7+1 ]; 
                                       
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1));
    

%% Shape description for N
material.parameters{3}.name                 = 'Ni';
material.parameters{3}.color_map            = [0 232 0]/255;
material.parameters{3}.sphere_centers       = [             0.5 ,            0.5 ,           0.5 ;
                                                  1-material.x5 ,    material.x5 ,   material.x5 ;
                                                    material.x5 ,  1-material.x5 ,   material.x5 ;
                                                    material.x5 ,    material.x5 , 1-material.x5 ; 
                                                    material.x5 ,  1-material.x5 , 1-material.x5 ;
                                                  1-material.x5 ,    material.x5 , 1-material.x5 ;
                                                  1-material.x5 ,  1-material.x5 ,   material.x5 ];
                                                                                          
material.parameters{3}.sphere_radius        = sphere_radius(3)*ones(1,size(material.parameters{3}.sphere_centers,1));
    

%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map, material.parameters{3}.color_map};

    