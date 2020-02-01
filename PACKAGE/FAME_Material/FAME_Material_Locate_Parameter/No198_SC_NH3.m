function material = No198_SC_NH3(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No198_SC_NH3
%
% See detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A3B_cP16_198_b_a.html
%
% Edit at 2017/6/20 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'NH3';
material.material_num = 2;
material.lattice_type = 'simple_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A3B_cP16_198_b_a.html';
material.x1           = 0.2107 ;
material.x2           = 0.3689 ;
material.y2           = 0.2671 ;
material.z2           = 0.1159 ;

material.lattice_constant.a = 5.13050;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for N
material.parameters{1}.name                 = 'N';
material.parameters{1}.color_map            = [0 0 255]/255;
material.parameters{1}.sphere_centers       = [        material.x1 ,      material.x1 ,        material.x1 ;
                                                   0.5-material.x1 ,    1-material.x1 ,    0.5+material.x1 ;
                                                     1-material.x1 ,  0.5+material.x1 ,    0.5-material.x1 ;
                                                   0.5+material.x1 ,  0.5-material.x1 ,      1-material.x1 ];
                         
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1));

%% Shape description for H
material.parameters{2}.name                 = 'H';
material.parameters{2}.color_map            = [255 255 255]/255;
material.parameters{2}.sphere_centers       = [       material.x2 ,      material.y2 ,        material.z2 ; % 5
                                                  0.5-material.x2 ,    1-material.y2 ,    0.5+material.z2 ; % 6
                                                    1-material.x2 ,  0.5+material.y2 ,    0.5-material.z2 ; % 7
                                                  0.5+material.x2 ,  0.5-material.y2 ,      1-material.z2 ; % 8
                                                      material.z2 ,      material.x2 ,        material.y2 ; % 9
                                                  0.5-material.z2 ,    1-material.x2 ,    0.5+material.y2 ; % 10
                                                    1-material.z2 ,  0.5+material.x2 ,    0.5-material.y2 ; % 11
                                                  0.5+material.z2 ,  0.5-material.x2 ,      1-material.y2 ; % 12
                                                      material.y2 ,      material.z2 ,        material.x2 ; % 13
                                                  0.5-material.y2 ,    1-material.z2 ,    0.5+material.x2 ; % 14
                                                    1-material.y2 ,  0.5+material.z2 ,    0.5-material.x2 ; % 15
                                                  0.5+material.y2 ,  0.5-material.z2 ,      1-material.x2 ]; % 16
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1));                                              

%%
material.parameters{1}.cylinder_bot_centers =  [       material.x1 ,      material.x1 ,        material.x1 ;
                                                       material.x1 ,      material.x1 ,        material.x1 ;
                                                       material.x1 ,      material.x1 ,        material.x1 ;
                                                   0.5-material.x1 ,    1-material.x1 ,    0.5+material.x1 ;
                                                   0.5-material.x1 ,    1-material.x1 ,    0.5+material.x1 ;
                                                   0.5-material.x1 ,    1-material.x1 ,    0.5+material.x1 ;
                                                     1-material.x1 ,  0.5+material.x1 ,    0.5-material.x1 ;
                                                     1-material.x1 ,  0.5+material.x1 ,    0.5-material.x1 ;
                                                     1-material.x1 ,  0.5+material.x1 ,    0.5-material.x1 ;
                                                   0.5+material.x1 ,  0.5-material.x1 ,       1-material.x1;
                                                   0.5+material.x1 ,  0.5-material.x1 ,       1-material.x1;
                                                   0.5+material.x1 ,  0.5-material.x1 ,       1-material.x1];
material.parameters{1}.cylinder_top_centers = 0.5*[    material.x1+material.x2 ,      material.x1+material.y2 ,        material.x1+material.z2 ;
                                                       material.x1+material.z2 ,      material.x1+material.x2 ,        material.x1+material.y2 ;
                                                       material.x1+material.y2 ,      material.x1+material.z2 ,        material.x1+material.x2 ;
                                                   0.5-material.x1+0.5-material.x2 ,    1-material.x1+1-material.y2 ,    0.5+material.x1+ 0.5+material.z2 ;
                                                   0.5-material.x1+0.5-material.z2 ,    1-material.x1+1-material.x2 ,    0.5+material.x1+0.5+material.y2 ;
                                                   0.5-material.x1+0.5-material.y2 ,    1-material.x1+1-material.z2 ,    0.5+material.x1+0.5+material.x2 ;
                                                     1-material.x1+1-material.x2 ,  0.5+material.x1+0.5+material.y2 ,    0.5-material.x1+ 0.5-material.z2 ;
                                                     1-material.x1+ 1-material.z2 ,  0.5+material.x1+0.5+material.x2 ,    0.5-material.x1+0.5-material.y2 ;
                                                     1-material.x1+ 1-material.y2 ,  0.5+material.x1+0.5+material.z2 ,    0.5-material.x1+ 0.5-material.x2 ;
                                                   0.5+material.x1+0.5+material.x2 ,  0.5-material.x1+0.5-material.y2 ,       1-material.x1+1-material.z2;
                                                   0.5+material.x1+0.5+material.z2 ,  0.5-material.x1+0.5-material.x2 ,       1-material.x1+1-material.y2;
                                                   0.5+material.x1+ 0.5+material.y2 ,  0.5-material.x1+0.5-material.z2 ,       1-material.x1+1-material.x2     ];
material.parameters{1}.cylinder_radius        = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_top_centers,1));   
material.parameters{2}.cylinder_bot_centers = [ material.x2 ,      material.y2 ,        material.z2 ;
                                                material.z2 ,      material.x2 ,        material.y2 ;
                                                material.y2 ,      material.z2 ,        material.x2 ;
                                                0.5-material.x2 ,    1-material.y2 ,    0.5+material.z2 ;
                                                 0.5-material.z2 ,    1-material.x2 ,    0.5+material.y2 ;
                                                 0.5-material.y2 ,    1-material.z2 ,    0.5+material.x2 ;
                                                  1-material.x2 ,  0.5+material.y2 ,    0.5-material.z2 ;
                                                  1-material.z2 ,  0.5+material.x2 ,    0.5-material.y2 ;
                                                   1-material.y2 ,  0.5+material.z2 ,    0.5-material.x2 ;
                                                  0.5+material.x2 ,  0.5-material.y2 ,      1-material.z2 ;
                                                    0.5+material.z2 ,  0.5-material.x2 ,      1-material.y2 ;
                                                   0.5+material.y2 ,  0.5-material.z2 ,      1-material.x2 ];  
material.parameters{2}.cylinder_top_centers = material.parameters{1}.cylinder_top_centers;

material.parameters{2}.cylinder_radius        = cylinder_radius(2)*ones(1,size(material.parameters{2}.cylinder_top_centers,1));

%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};

    