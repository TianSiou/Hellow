function material = No012_Primitive_Monoclinic_beta_Ga2O3(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No012_Primitive_Monoclinic_beta_Ga2O3
%
% More detail on following website:
%        https://materialsproject.org/materials/mp-886/#
%
% Edit at 2017/8/16 By Hsiao-Han Huang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
material.Prototype    = 'beta_Ga2O3';
material.material_num = 2;

material.lattice_type = 'primitive_monoclinic';
material.website      = 'https://materialsproject.org/materials/mp-886/#';

material.lattice_constant.a = 3.08296874;
material.lattice_constant.b = 5.87615337;
material.lattice_constant.c = 12.45244640;
material.lattice_constant.alpha = pi*(180-103.68360198)/180;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );

%% Shape description for Ga (Sphere)
material.parameters{1}.name                 = 'Ga';
material.parameters{1}.color_map            = [128 255 255]/255;
material.parameters{1}.sphere_centers       = [   0.658410  0.000000  0.314082  ;
                                                  0.841590  0.500000  0.685918  ;
                                                  0.589878  0.500000  0.794761  ;
                                                  0.910122  0.000000  0.205239  ;
                                                  0.158410  0.500000  0.314082  ;
                                                  0.341590  0.000000  0.685918  ;
                                                  0.089878  0.000000  0.794761  ;
                                                  0.410122  0.500000  0.205239  ];
material.parameters{1}.sphere_centers(:,1)   = 1 - material.parameters{1}.sphere_centers(:,1);
material.parameters{1}.sphere_centers(:,2:3) = material.parameters{1}.sphere_centers(:,[3,2]);
% material.parameters{1}.sphere_centers = material.parameters{1}.sphere_centers(:,[2,1,3]);                                                   
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(size(material.parameters{1}.sphere_centers,1),1);                                              

%% Shape description for O (Sphere)
material.parameters{2}.name                 = 'O';
material.parameters{2}.color_map            = [255 0 0]/255;

material.parameters{2}.sphere_centers   = [       0.995897  0.500000  0.256491  ;
                                                  0.504103  0.000000  0.743509  ;
                                                  0.673598  0.500000  0.564294  ;
                                                  0.826402  0.000000  0.435706  ;
                                                  0.836497  0.000000  0.891047  ;
                                                  0.663503  0.500000  0.108953  ;
                                                  0.495897  0.000000  0.256491  ;
                                                  0.004103  0.500000  0.743509  ;
                                                  0.173598  0.000000  0.564294  ;
                                                  0.326402  0.500000  0.435706  ;
                                                  0.336497  0.500000  0.891047  ;
                                                  0.163503  0.000000  0.108953  ];
material.parameters{2}.sphere_centers(:,1)   = 1 - material.parameters{2}.sphere_centers(:,1);
material.parameters{2}.sphere_centers(:,2:3) = material.parameters{2}.sphere_centers(:,[3,2]);                                 
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(size(material.parameters{2}.sphere_centers,1),1);                                              

%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};