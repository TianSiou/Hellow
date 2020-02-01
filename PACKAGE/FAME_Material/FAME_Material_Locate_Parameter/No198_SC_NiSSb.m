function material = No198_SC_NiSSb(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No198_SC_NiSSb
%
% Other compounds with this structure:AsBaPt, AsPdS, BiIrS, BiRhSe, CaPtSi, CrPtSb, EuPtSi, IrLaSi, IrSbSe, many others
% See detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/ABC_cP12_198_a_a_a.html
% Edit at 2017/6/20 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
material.Prototype    = 'NiSSb';
material.material_num = 3;

material.lattice_type = 'simple_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/ABC_cP12_198_a_a_a.html';
material.x1           = -0.024 ;
material.x2           =   0.39 ;
material.x3           =  0.875 ;

material.lattice_constant.a = 5.881;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Ni
material.parameters{1}.name                 = 'Ni';
material.parameters{1}.color_map            = [76 199 76]/255;
material.parameters{1}.sphere_centers       = [     1+material.x1 ,   1+material.x1 ,    1+material.x1 ;  % B1
                                                  0.5-material.x1 ,    -material.x1 ,  0.5+material.x1 ;  % B2
                                                     -material.x1 , 0.5+material.x1 ,  0.5-material.x1 ;  % B3
                                                  0.5+material.x1 , 0.5-material.x1 ,     -material.x1 ]; % B4

material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 
material.parameters{1}.cylinder_bot_centers =  [   0.5-material.x1 ,    -material.x1 ,  0.5+material.x1 ;
                                                   0.5-material.x1 ,    -material.x1 ,  0.5+material.x1 ;
                                                      -material.x1 , 0.5+material.x1 ,  0.5-material.x1 ;
                                                      -material.x1 , 0.5+material.x1 ,  0.5-material.x1 ;
                                                   0.5+material.x1 , 0.5-material.x1 ,      -material.x1;
                                                   0.5+material.x1 , 0.5-material.x1 ,      -material.x1];
                                               
material.parameters{1}.cylinder_top_centers = 0.5*[     0.5-material.x1+material.x2 ,        -material.x1+material.x2 ,      0.5+material.x1+material.x2 ;  % 2&5
                                                    0.5-material.x1+0.5+material.x2 ,    -material.x1+0.5-material.x2 ,    0.5+material.x1-material.x2+1 ;  % 2&8
                                                           -material.x1+material.x2 ,     0.5+material.x1+material.x2 ,      0.5-material.x1+material.x2 ;  % 3&5
                                                       -material.x1+0.5-material.x2 ,   0.5+material.x1-material.x2+1 ,  0.5-material.x1+0.5+material.x2 ;  % 3&6                             
                                                        0.5+material.x1+material.x2 ,     0.5-material.x1+material.x2 ,          -material.x1+material.x2;  % 4&5
                                                      0.5+material.x1-material.x2+1 , 0.5-material.x1+0.5+material.x2 ,      -material.x1+0.5-material.x2]; % 4&7
                                  

material.parameters{1}.cylinder_radius        = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_bot_centers,1)); 
%% Shape description for S
material.parameters{2}.name                 = 'S';
material.parameters{2}.color_map            = [246 246 44]/255;
material.parameters{2}.sphere_centers       = [        material.x2 ,      material.x2 ,        material.x2 ;  % B5
                                                   0.5-material.x2 ,   -material.x2+1 ,    0.5+material.x2 ;  % B6
                                                    -material.x2+1 ,  0.5+material.x2 ,    0.5-material.x2 ;  % B7
                                                   0.5+material.x2 ,  0.5-material.x2 ,     -material.x2+1 ]; % B8

material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 
 material.parameters{2}.cylinder_bot_centers =  [      material.x2 ,      material.x2 ,        material.x2 ;
                                                       material.x2 ,      material.x2 ,        material.x2 ;
                                                       material.x2 ,      material.x2 ,        material.x2 ;
                                                       material.x2 ,      material.x2 ,        material.x2 ;
                                                       material.x2 ,      material.x2 ,        material.x2 ;
                                                       material.x2 ,      material.x2 ,        material.x2 ;
                                                   0.5-material.x2 ,   -material.x2+1 ,    0.5+material.x2 ;
                                                   0.5-material.x2 ,   -material.x2+1 ,    0.5+material.x2 ;
                                                    -material.x2+1 ,  0.5+material.x2 ,    0.5-material.x2 ;
                                                    -material.x2+1 ,  0.5+material.x2 ,    0.5-material.x2 ;
                                                   0.5+material.x2 ,  0.5-material.x2 ,      -material.x2+1;
                                                   0.5+material.x2 ,  0.5-material.x2 ,      -material.x2+1];
                                                   
material.parameters{2}.cylinder_top_centers = 0.5*[    material.x2+0.5-material.x1 ,          material.x2-material.x1 ,       material.x2+0.5+material.x1 ;   % 5&2
                                                           material.x2-material.x1 ,      material.x2+0.5+material.x1 ,       material.x2+0.5-material.x1 ;   % 5&3
                                                       material.x2+0.5+material.x1 ,      material.x2+0.5-material.x1 ,           material.x2-material.x1 ;   % 5&4
                                                     material.x2+0.5-material.x3+1 ,        material.x2-material.x3+1 ,     material.x2+0.5+material.x3-1 ;   % 5&10
                                                         material.x2-material.x3+1 ,    material.x2+0.5+material.x3-1 ,     material.x2+0.5-material.x3+1 ;   % 5&11
                                                     material.x2+0.5+material.x3-1 ,    material.x2+0.5-material.x3+1 ,         material.x2-material.x3+1 ;   % 5&12
                                                      -material.x1+0.5-material.x2 ,    0.5+material.x1-material.x2+1 ,   0.5-material.x1+0.5+material.x2 ;   % 6&3
                                                    -material.x3+1+0.5-material.x2 ,  0.5+material.x3-1-material.x2+1 , 0.5-material.x3+1+0.5+material.x2 ;   % 6&11
                                                     0.5+material.x1-material.x2+1 ,  0.5-material.x1+0.5+material.x2 ,       -material.x1+0.5-material.x2;   % 7&4
                                                   0.5+material.x3-1-material.x2+1 ,0.5-material.x3+1+0.5+material.x2 ,    -material.x3+1+0.5-material.x2 ;   % 7&12
                                                   0.5-material.x1+0.5+material.x2 ,     -material.x1+0.5-material.x2 ,    0.5+material.x1-material.x2+1  ;   % 8&2
                                                 0.5-material.x3+1+0.5+material.x2 ,   -material.x3+1+0.5-material.x2 ,  0.5+material.x3-1-material.x2+1  ];  % 8&10

material.parameters{2}.cylinder_radius = cylinder_radius(2)*ones(1,size(material.parameters{2}.cylinder_bot_centers,1)); 
%% Shape description for Sb
material.parameters{3}.name                 = 'Sb';
material.parameters{3}.color_map            = [150 93 172]/255;
material.parameters{3}.sphere_centers       = [         material.x3 ,        material.x3 ,          material.x3 ;  % B9
                                                  0.5-material.x3+1 ,     -material.x3+1 ,    0.5+material.x3-1 ;  % B10
                                                     -material.x3+1 ,  0.5+material.x3-1 ,    0.5-material.x3+1 ;  % B11
                                                  0.5+material.x3-1 ,  0.5-material.x3+1 ,       -material.x3+1 ]; % B12

material.parameters{3}.sphere_radius        = sphere_radius(3)*ones(1,size(material.parameters{3}.sphere_centers,1)); 

material.parameters{3}.cylinder_bot_centers =  [   0.5-material.x3+1 ,     -material.x3+1 ,    0.5+material.x3-1 ;
                                                   0.5-material.x3+1 ,     -material.x3+1 ,    0.5+material.x3-1 ;
                                                      -material.x3+1 ,  0.5+material.x3-1 ,    0.5-material.x3+1 ;
                                                      -material.x3+1 ,  0.5+material.x3-1 ,    0.5-material.x3+1 ;
                                                   0.5+material.x3-1 ,  0.5-material.x3+1 ,       -material.x3+1 ;
                                                   0.5+material.x3-1 ,  0.5-material.x3+1 ,       -material.x3+1 ];
                                               
material.parameters{3}.cylinder_top_centers = 0.5*[     0.5-material.x3+1+material.x2 ,        -material.x3+1+material.x2 ,       0.5+material.x3-1+material.x2 ;   % 10&5
                                                    0.5-material.x3+1+0.5+material.x2 ,    -material.x3+1+0.5-material.x2 ,     0.5+material.x3-1-material.x2+1 ;   % 10&8
                                                           -material.x3+1+material.x2 ,     0.5+material.x3-1+material.x2 ,       0.5-material.x3+1+material.x2 ;   % 11&5
                                                       -material.x3+1+0.5-material.x2 ,   0.5+material.x3-1-material.x2+1 ,   0.5-material.x3+1+0.5+material.x2 ;   % 11&6
                                                        0.5+material.x3-1+material.x2 ,     0.5-material.x3+1+material.x2 ,          -material.x3+1+material.x2 ;   % 12&5
                                                      0.5+material.x3-1-material.x2+1 , 0.5-material.x3+1+0.5+material.x2 ,      -material.x3+1+0.5-material.x2 ];  % 12&7
                                                  
material.parameters{3}.cylinder_radius = cylinder_radius(3)*ones(1,size(material.parameters{3}.cylinder_bot_centers,1)); 
%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map,material.parameters{3}.color_map};