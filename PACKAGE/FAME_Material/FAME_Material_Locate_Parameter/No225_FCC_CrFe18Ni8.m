function material = No225_FCC_CrFe18Ni8(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No225_FCC_CrFe18Ni8
%
% Austenitic steels are alloys of iron and other metals with an averaged face-centered cubic structure. 
% If we set x2=1/3x2=1/3, x3=2/3x3=2/3, and y4=2/3y4=2/3, the atoms are on the sites of an fcc lattice with lattice constant afcc=1/3a.
%
% See detail on following website:
%      http://www.aflowlib.org/CrystalDatabase/AB18C8_cF108_225_a_eh_f.html
% Edit at 2017/6/22 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'CrFe18Ni8';
material.material_num = 3;

material.lattice_type = 'face_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB18C8_cF108_225_a_eh_f.html';
material.x2           = 0.325 ;
material.x3           = 0.65833 ;
material.y4           = 0.66 ;
material.lattice_constant.a = 10.56000;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Cr
material.parameters{1}.name                 = 'Cr';
material.parameters{1}.color_map            = [104 104 255]/255;
material.parameters{1}.sphere_centers       = [0, 0, 0];
                                          
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1));  

%% Shape description for Fe
material.parameters{2}.name                 = 'Fe';
material.parameters{2}.color_map            =  [188 96 69]/255;
material.parameters{2}.sphere_centers       = [       -material.x2+1 ,        material.x2 ,        material.x2 ;
                                                         material.x2 ,      1-material.x2 ,        material.x2 ;
                                                         material.x2 ,        material.x2 ,      1-material.x2 ;
                                                         material.x2 ,     -material.x2+1 ,     -material.x2+1 ;
                                                      -material.x2+1 ,        material.x2 ,     -material.x2+1 ;
                                                      -material.x2+1 ,     -material.x2+1 ,        material.x2 ;
                                                     2*material.y4-1 ,                  0 ,                  0 ;
                                                                   0 ,    2*material.y4-1 ,   -2*material.y4+2 ;
                                                                   0 ,   -2*material.y4+2 ,    2*material.y4-1 ;
                                                    -2*material.y4+2 ,                  0 ,                  0 ;
                                                                   0 ,    2*material.y4-1 ,                  0 ;
                                                    -2*material.y4+2 ,                  0 ,    2*material.y4-1 ;
                                                     2*material.y4-1 ,                  0 ,   -2*material.y4+2 ;
                                                                   0 ,   -2*material.y4+2 ,                  0 ;
                                                                   0 ,                  0 ,    2*material.y4-1 ;
                                                     2*material.y4-1 ,   -2*material.y4+2 ,                  0 ;
                                                    -2*material.y4+2 ,    2*material.y4-1 ,                  0 ;
                                                                   0 ,                  0 ,   -2*material.y4+2 ];
                                                                                              
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1));  
   
%% Shape description for Ni
material.parameters{3}.name                 = 'Ni';
material.parameters{3}.color_map            = [0 232 0]/255;
material.parameters{3}.sphere_centers       =[         material.x3 ,       material.x3 ,       material.x3 ;
                                                       material.x3 ,       material.x3 ,  -3*material.x3+2 ;
                                                       material.x3 , -3* material.x3+2 ,       material.x3 ;
                                                 -3* material.x3+2 ,       material.x3 ,       material.x3 ;
                                                    -material.x3+1 ,    -material.x3+1 ,   3*material.x3-1 ;
                                                    -material.x3+1 ,    -material.x3+1 ,    -material.x3+1 ;
                                                    -material.x3+1 ,  3* material.x3-1 ,    -material.x3+1 ;
                                                  3* material.x3-1 ,    -material.x3+1 ,    -material.x3+1 ];                                         
material.parameters{3}.sphere_radius        = sphere_radius(3)*ones(1,size(material.parameters{3}.sphere_centers,1));  
   
%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map, material.parameters{3}.color_map};

    