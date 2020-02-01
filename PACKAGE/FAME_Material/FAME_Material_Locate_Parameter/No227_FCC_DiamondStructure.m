function material = No227_FCC_DiamondStructure(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No227_FCC_DiamondStructure
%
% Edit at 2017/6/25 By Jia-Wei Lin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'Si';
material.material_num = 1;

material.lattice_type = 'face_centered_cubic';
material.website      = [];

material.lattice_constant.a = 1;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Si (Sphere)
material.parameters{1}.name                 = 'Si';
material.parameters{1}.color_map            = [138 128 119]/255;
material.parameters{1}.sphere_centers       = [    0,    0,    0;
                                                0.25, 0.25, 0.25 ];                                             
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1));

                                      
%% %% Shape description for Si (Cylinder)
material.parameters{1}.cylinder_bot_centers =  [ 0,  0,  0;
                                                 1,  0,  0;
                                                 0,  1,  0;
                                                 0,  0,  1];
material.parameters{1}.cylinder_top_centers = [  material.parameters{1}.sphere_centers(2,:);
                                                 material.parameters{1}.sphere_centers(2,:);
                                                 material.parameters{1}.sphere_centers(2,:);
                                                 material.parameters{1}.sphere_centers(2,:)];
                                             
material.parameters{1}.cylinder_radius        = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_top_centers,1));

%% Collect color map
material.color_map = {material.parameters{1}.color_map};