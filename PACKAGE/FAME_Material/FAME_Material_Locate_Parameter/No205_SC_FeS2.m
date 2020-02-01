function material = No205_SC_FeS2(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No205_SC_FeS2
%
% Gives crystalline data for weakly anisotropic pyrite which we have tabulated as P1 FeS2. Also gives crystallographic data for the cubic pyrite structure. 
% Also see the C18 (marcasite) FeS2 structure.
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/AB2_cP12_205_a_c.html
%
% Edit at 2017/7/7 By Hsiao-Han Huang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'FeS2';
material.material_num = 2;

material.lattice_type = 'simple_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB2_cP12_205_a_c.html';
material.x2           = 0.3851;

material.lattice_constant.a = 5.417;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Fe
material.parameters{1}.name                 = 'Fe';
material.parameters{1}.color_map            = [124 6 3]/255;
material.parameters{1}.sphere_centers       = [     0,     0,    0; 
                                                   .5,     0,   .5;
                                                    0,    .5,   .5;
                                                   .5,    .5,    0];
                                               
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 

material.parameters{1}.cylinder_bot_centers = [     .5,      0,     .5;
                                                     0,     .5,     .5;
                                                    .5,     .5,      0;
                                                     0,     .5,     .5;
                                                    .5,     .5,      1;
                                                     0,      1,      1;
                                                    .5,     .5,      0;
                                                     1,      1,      0;
                                                    .5,      1,     .5;
                                                    .5,      0,     .5;
                                                     1,     .5,     .5;
                                                     1,      0,      1;
                                                    .5,     .5,      1;
                                                    .5,      1,     .5;
                                                     1,     .5,     .5;
                                                    .5,     .5,      0;
                                                     1,      0,      0;
                                                     1,     .5,     .5;
                                                     0,      0,      1;
                                                    .5,     .5,      1;
                                                    .5,      0,     .5;
                                                     0,      1,      0;
                                                     0,     .5,     .5;
                                                    .5,      1,     .5]; 
material.parameters{1}.cylinder_top_centers = [     .5*(material.x2+.5),           .5*material.x2,       .5*(material.x2+.5);
                                                         .5*material.x2,      .5*(material.x2+.5),       .5*(material.x2+.5);
                                                    .5*(material.x2+.5),      .5*(material.x2+.5),            .5*material.x2;
                                                   .5*(-material.x2+.5),    .5*(-material.x2+1.5),        .5*(material.x2+1);
                                                    .5*(-material.x2+1),    .5*(-material.x2+1.5),      .5*(material.x2+1.5);
                                                   .5*(-material.x2+.5),      .5*(-material.x2+2),      .5*(material.x2+1.5);
                                                  .5*(-material.x2+1.5),       .5*(material.x2+1),      .5*(-material.x2+.5);
                                                    .5*(-material.x2+2),     .5*(material.x2+1.5),      .5*(-material.x2+.5);
                                                  .5*(-material.x2+1.5),     .5*(material.x2+1.5),       .5*(-material.x2+1);
                                                     .5*(material.x2+1),     .5*(-material.x2+.5),     .5*(-material.x2+1.5);
                                                   .5*(material.x2+1.5),      .5*(-material.x2+1),     .5*(-material.x2+1.5);
                                                   .5*(material.x2+1.5),     .5*(-material.x2+.5),       .5*(-material.x2+2);
                                                  .5*(-material.x2+1.5),    .5*(-material.x2+1.5),       .5*(-material.x2+2);
                                                  .5*(-material.x2+1.5),      .5*(-material.x2+2),     .5*(-material.x2+1.5);
                                                    .5*(-material.x2+2),    .5*(-material.x2+1.5),     .5*(-material.x2+1.5);
                                                     .5*(material.x2+1),      .5*(material.x2+.5),      .5*(-material.x2+.5);
                                                   .5*(material.x2+1.5),         .5*(material.x2),      .5*(-material.x2+.5);
                                                   .5*(material.x2+1.5),      .5*(material.x2+.5),       .5*(-material.x2+1);
                                                       .5*(material.x2),     .5*(-material.x2+.5),      .5*(material.x2+1.5);
                                                    .5*(material.x2+.5),      .5*(-material.x2+1),      .5*(material.x2+1.5);
                                                    .5*(material.x2+.5),     .5*(-material.x2+.5),        .5*(material.x2+1);
                                                   .5*(-material.x2+.5),     .5*(material.x2+1.5),          .5*(material.x2);
                                                   .5*(-material.x2+.5),       .5*(material.x2+1),       .5*(material.x2+.5);
                                                    .5*(-material.x2+1),     .5*(material.x2+1.5),       .5*(material.x2+.5);
                                                   ];
                                           
material.parameters{1}.cylinder_radius        = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_bot_centers,1)); 
%% Shape description for S
material.parameters{2}.name                 = 'S';
material.parameters{2}.color_map            = [214 218 37]/255;
material.parameters{2}.sphere_centers       = [      material.x2,       material.x2,       material.x2;
                                                  .5-material.x2,    -material.x2+1,    .5+material.x2;
                                                  -material.x2+1,    .5+material.x2,    .5-material.x2;
                                                  .5+material.x2,    .5-material.x2,    -material.x2+1;
                                                  -material.x2+1,    -material.x2+1,    -material.x2+1;
                                                  .5+material.x2,       material.x2,    .5-material.x2;
                                                     material.x2,    .5-material.x2,    .5+material.x2;
                                                  .5-material.x2,    .5+material.x2,       material.x2];
                                           
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 

material.parameters{2}.cylinder_bot_centers = [     material.x2,       material.x2,       material.x2;
                                                    material.x2,       material.x2,       material.x2;
                                                    material.x2,       material.x2,       material.x2;
                                                 .5-material.x2,    -material.x2+1,    .5+material.x2;
                                                 .5-material.x2,    -material.x2+1,    .5+material.x2;
                                                 .5-material.x2,    -material.x2+1,    .5+material.x2;
                                                 -material.x2+1,    .5+material.x2,    .5-material.x2;
                                                 -material.x2+1,    .5+material.x2,    .5-material.x2;
                                                 -material.x2+1,    .5+material.x2,    .5-material.x2;
                                                 .5+material.x2,    .5-material.x2,    -material.x2+1;
                                                 .5+material.x2,    .5-material.x2,    -material.x2+1;
                                                 .5+material.x2,    .5-material.x2,    -material.x2+1;
                                                 -material.x2+1,    -material.x2+1,    -material.x2+1;
                                                 -material.x2+1,    -material.x2+1,    -material.x2+1;
                                                 -material.x2+1,    -material.x2+1,    -material.x2+1;
                                                 .5+material.x2,       material.x2,    .5-material.x2;
                                                 .5+material.x2,       material.x2,    .5-material.x2;
                                                 .5+material.x2,       material.x2,    .5-material.x2;
                                                    material.x2,    .5-material.x2,    .5+material.x2;
                                                    material.x2,    .5-material.x2,    .5+material.x2;
                                                    material.x2,    .5-material.x2,    .5+material.x2;
                                                 .5-material.x2,    .5+material.x2,       material.x2;
                                                 .5-material.x2,    .5+material.x2,       material.x2;
                                                 .5-material.x2,    .5+material.x2,       material.x2;
                                                      ];
                                              
material.parameters{2}.cylinder_top_centers = [    .5*(material.x2+.5),           .5*material.x2,      .5*(material.x2+.5);
                                                        .5*material.x2,      .5*(material.x2+.5),      .5*(material.x2+.5);
                                                   .5*(material.x2+.5),      .5*(material.x2+.5),           .5*material.x2;
                                                  .5*(-material.x2+.5),    .5*(-material.x2+1.5),       .5*(material.x2+1);
                                                   .5*(-material.x2+1),    .5*(-material.x2+1.5),     .5*(material.x2+1.5);
                                                  .5*(-material.x2+.5),      .5*(-material.x2+2),     .5*(material.x2+1.5);
                                                 .5*(-material.x2+1.5),       .5*(material.x2+1),     .5*(-material.x2+.5);
                                                   .5*(-material.x2+2),     .5*(material.x2+1.5),     .5*(-material.x2+.5);
                                                 .5*(-material.x2+1.5),     .5*(material.x2+1.5),      .5*(-material.x2+1);
                                                    .5*(material.x2+1),     .5*(-material.x2+.5),    .5*(-material.x2+1.5);
                                                  .5*(material.x2+1.5),      .5*(-material.x2+1),    .5*(-material.x2+1.5);
                                                  .5*(material.x2+1.5),     .5*(-material.x2+.5),      .5*(-material.x2+2);
                                                 .5*(-material.x2+1.5),    .5*(-material.x2+1.5),      .5*(-material.x2+2);
                                                 .5*(-material.x2+1.5),      .5*(-material.x2+2),    .5*(-material.x2+1.5);
                                                   .5*(-material.x2+2),    .5*(-material.x2+1.5),    .5*(-material.x2+1.5);
                                                    .5*(material.x2+1),      .5*(material.x2+.5),     .5*(-material.x2+.5);
                                                  .5*(material.x2+1.5),         .5*(material.x2),     .5*(-material.x2+.5);
                                                  .5*(material.x2+1.5),      .5*(material.x2+.5),      .5*(-material.x2+1);
                                                      .5*(material.x2),     .5*(-material.x2+.5),     .5*(material.x2+1.5);
                                                   .5*(material.x2+.5),      .5*(-material.x2+1),     .5*(material.x2+1.5);
                                                   .5*(material.x2+.5),     .5*(-material.x2+.5),       .5*(material.x2+1);
                                                  .5*(-material.x2+.5),     .5*(material.x2+1.5),         .5*(material.x2);
                                                  .5*(-material.x2+.5),       .5*(material.x2+1),      .5*(material.x2+.5);
                                                   .5*(-material.x2+1),     .5*(material.x2+1.5),      .5*(material.x2+.5);
                                                   ];
                                           
material.parameters{2}.cylinder_radius        = cylinder_radius(2)*ones(1,size(material.parameters{2}.cylinder_bot_centers,1)); 
%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};
