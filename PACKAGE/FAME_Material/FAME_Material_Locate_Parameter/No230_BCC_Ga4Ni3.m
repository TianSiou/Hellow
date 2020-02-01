function material = No230_BCC_Ga4Ni3(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No230_BCC_Ga4Ni3
% This is a simple defect superstructure of the CsCl (B2) structure. If a GaNi B2 
%structure is expanded into a 128 atom supercell, we can describe it using space group
%Ia3d (#230), with Ga atoms on the (16a) and (48f) Wyckoff sites and Ni atoms on the (16b)
%and (48g) sites. Removing the Ni atoms from the (16b) sites yields this structure.
%
% See detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A4B3_cI112_230_af_g.html
% Edit at 2017/6/22 By Jyun-Wei Lin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'Ga4Ni3';
material.material_num = 2;

material.lattice_type = 'body_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/A4B3_cI112_230_af_g.html';
material.lattice_constant.a           = 11.411  ;
material.x2          = 0       ;
material.y3          = 0.625   ;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Ga
material.parameters{1}.name           = 'Ga';
material.parameters{1}.color_map      = [172 125 200]/255;
material.parameters{1}.sphere_centers = [               0 ,                  0 ,                 0 ;
                                                      0.5 ,                  0 ,               0.5 ;
                                                        0 ,                0.5 ,               0.5 ;
                                                      0.5 ,                0.5 ,                 0 ;
                                                      0.5 ,                  0 ,                 0 ;
                                                      0.5 ,                0.5 ,               0.5 ;
                                                        0 ,                  0 ,               0.5 ;
                                                        0 ,                0.5 ,                 0 ;
                                                     0.25 ,   0.25+material.x2 ,       material.x2 ;
                                                     0.75 ,   0.25-material.x2 ,   0.5-material.x2 ;
                                              material.x2 ,               0.25 ,  0.25+material.x2 ;
                                          0.5-material.x2 ,               0.75 ,  0.25-material.x2 ;
                                         0.25+material.x2 ,        material.x2 ,              0.25 ;
                                         0.25-material.x2 ,    0.5-material.x2 ,              0.75 ;
                                         0.25+material.x2 ,               0.75 ,       material.x2 ;
                                         0.25-material.x2 ,               0.25 ,   0.5-material.x2 ;
                                                     0.75 ,        material.x2 ,  0.25+material.x2 ;
                                                     0.25 ,    0.5-material.x2 ,  0.25-material.x2 ;
                                          0.5-material.x2 ,   0.25-material.x2 ,              0.25 ;
                                              material.x2 ,   0.25+material.x2 ,              0.75 ;
                                                     0.75 ,   0.75-material.x2 ,      -material.x2 ;
                                                     0.25 ,   0.75+material.x2 ,   0.5+material.x2 ;
                                             -material.x2 ,               0.75 ,  0.75-material.x2 ;
                                          0.5+material.x2 ,               0.25 ,  0.75+material.x2 ;
                                         0.75-material.x2 ,       -material.x2 ,              0.75 ;
                                         0.75+material.x2 ,    0.5+material.x2 ,              0.25 ;
                                         0.75-material.x2 ,               0.25 ,      -material.x2 ;  
                                         0.75+material.x2 ,               0.75 ,   0.5+material.x2 ;
                                                     0.25 ,       -material.x2 ,  0.75-material.x2 ;
                                                     0.75 ,    0.5+material.x2 ,  0.75+material.x2 ;
                                          0.5+material.x2 ,   0.75+material.x2 ,              0.75 ;
                                             -material.x2 ,   0.75-material.x2 ,              0.25 ];
material.parameters{1}.sphere_radius= sphere_radius(1)*ones(size(material.parameters{1}.sphere_centers,1),1)';
%% Shape description for Ni
material.parameters{2}.name           = 'Ni';
material.parameters{2}.color_map      = [27 235 27]/255;
for i=1:4
    pt(2*i)=0.25*i+2*material.y3;
    nt(2*i)=0.25*i-2*material.y3;
    pt(2*i-1)=0.125*(2*i-1)+material.y3;
    nt(2*i-1)=0.125*(2*i-1)-material.y3;
end
material.parameters{2}.sphere_centers = [    0.25,  nt(3),  pt(1) ;
                                            nt(6),  nt(1),  nt(3) ;
                                            pt(6),  pt(1),  pt(3) ;
                                             0.25,  pt(3),  nt(1) ;
                                            pt(1),   0.25,  nt(3) ;
                                            nt(3),  nt(6),  nt(1) ;
                                            pt(3),  pt(6),  pt(1) ;
                                            nt(1),   0.25,  pt(3) ;
                                            nt(3),  pt(1),   0.25 ;
                                            nt(1),  nt(3),  nt(6) ;
                                            pt(1),  pt(3),  pt(6) ;
                                            pt(3),  nt(1),   0.25 ;
                                             0.75,  pt(5),  nt(7) ;
                                            pt(2),  pt(7),  pt(5) ;
                                            nt(2),  nt(7),  nt(5) ;
                                             0.75,  nt(5),  pt(7) ;
                                            nt(7),   0.75,  pt(5) ;
                                            pt(5),  pt(2),  pt(7) ;
                                            nt(5),  nt(2),  nt(7) ;
                                            pt(7),   0.75,  nt(5) ;
                                            pt(5),  nt(7),   0.75 ;
                                            pt(7),  pt(5),  pt(2) ;
                                            nt(7),  nt(5),  nt(2) ;
                                            nt(5),  pt(7),   0.75 ];
material.parameters{2}.sphere_centers=mod(material.parameters{2}.sphere_centers,1);
material.parameters{2}.sphere_radius= sphere_radius(2)*ones(size(material.parameters{2}.sphere_centers,1),1)';
%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};

    