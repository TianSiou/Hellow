function material = No225_FCC_AlCu2Mn(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No225_FCC_AlCu2Mn
%
% See detail on following website:
%      http://www.aflowlib.org/CrystalDatabase/AB2C_cF16_225_a_c_b.html
% Edit at 2017/6/22 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'AlCu2Mn';
material.material_num = 3;

material.lattice_type = 'face_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB2C_cF16_225_a_c_b.html';
material.lattice_constant.a = 5.95000;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Al
material.parameters{1}.name                 = 'Al';
material.parameters{1}.color_map            = [128 128 128]/255;
material.parameters{1}.sphere_centers       = [    0 ,  0 ,   0 ];
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1));
    

%% Shape description for Cu
material.parameters{2}.name                 = 'Cu';
material.parameters{2}.color_map            =  [185 122 87]/255;
material.parameters{2}.sphere_centers       = [    0.25 ,  0.25 ,  0.25 ;
                                                   0.75 ,  0.75 ,  0.75 ];
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1));
    

%% Shape description for Mn
material.parameters{3}.name                 = 'Mn';
material.parameters{3}.color_map            = [163 73 164]/255;
material.parameters{3}.sphere_centers       = [    0.5 ,  0.5 ,  0.5 ];
material.parameters{3}.sphere_radius        = sphere_radius(3)*ones(1,size(material.parameters{3}.sphere_centers,1));
    


%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map, material.parameters{3}.color_map};

    