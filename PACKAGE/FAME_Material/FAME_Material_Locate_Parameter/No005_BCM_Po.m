function material = No005_BCM_Po(sphere_radius,cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No005_BCM_Po
%
% This was the original determination of the structure of Po, and given the Strukturbericht 
% designation A19. (Gottfried, 1938, 4-5). Eventually it was determined that the sample 
% used here was a mixture of £\¡VPo (Ah) and £]¡VPo (Ai) (Donohue, 1982, 390). We retain 
% the A19 page for historical interest.
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A_mC12_5_3c.html
%
% Edit at 2017/7/5 By Jia-Wei Lin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'Po';
material.material_num = 1;

material.lattice_type = 'base_centered_monoclinic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A_mC12_5_3c.html';

material.x1           = 0.05;
material.y1           = 0.27;
material.z1           = 0.245;
material.x2           = 0.63;
material.y2           = 0.3;
material.z2           = 0.4;
material.x3           = 0.245;
material.y3           = 0.43;
material.z3           = 0.07;

material.lattice_constant.a       = 7.42;
material.lattice_constant.b       = 4.29;
material.lattice_constant.c       = 14.1;
material.lattice_constant.alpha   = pi*(180-92)/180;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );

%% Shape description for Po (Sphere)
material.parameters{1}.name                 = 'Po';
material.parameters{1}.color_map            = [135 35 188]/255;
% material.parameters{1}.sphere_centers       = [   1 + (material.x1-material.y1),      material.x1+material.y1,       material.z1        ; 
%                                                   1 - (material.x1+material.y1),      material.y1-material.x1,   1 - material.z1        ;
%                                                       (material.x2-material.y2),      material.x2+material.y2,       material.z2        ;
%                                                   1 - (material.x2+material.y2),  1 + material.y2-material.x2,   1 - material.z2        ;
%                                                   1 + (material.x3-material.y3),      material.x3+material.y3,       material.z3        ;
%                                                   1 - (material.x3+material.y3),      material.y3-material.x3,   1 - material.z3        ];
material.parameters{1}.sphere_centers       = [     material.x1+material.y1,    1 + (material.x1-material.y1),          material.z1        ; 
                                                    material.y1-material.x1,    1 - (material.x1+material.y1),      1 - material.z1        ;
                                                    material.x2+material.y2,        (material.x2-material.y2),          material.z2        ;
                                                1 + material.y2-material.x2,    1 - (material.x2+material.y2),      1 - material.z2        ;
                                                    material.x3+material.y3,    1 + (material.x3-material.y3),          material.z3        ;
                                                    material.y3-material.x3,    1 - (material.x3+material.y3),      1 - material.z3        ];
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1));

%% Collect color map
material.color_map = {material.parameters{1}.color_map};