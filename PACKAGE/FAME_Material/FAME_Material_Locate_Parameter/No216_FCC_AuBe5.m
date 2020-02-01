function material = No216_FCC_AuBe5(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No216_FCC_AuBe5
%
% The lattice constant for this structure is taken from (Batchelder, 1958),which does not give the internal coordinate for the (16c) site.
% However, (Baenziger, 1950) assumes that uranium compounds of this type have an internal parameter x3?5/8.(Pearson, 1958) uses this to infer a value of x3?5/8 here as well.
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/AB5_cF24_216_a_ce.html
%
% Edit at 2017/6/19 By Yu-Shiuan Jian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'AuBe5';
material.material_num = 2;

material.lattice_type = 'face_centered_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB5_cF24_216_a_ce.html';
material.x3           = 5/8;

material.lattice_constant.a = 6.1;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for C
material.parameters{1}.name                 = 'Au';
material.parameters{1}.color_map            = [243 222 35]/255;

material.parameters{1}.sphere_centers       = [ 0, 0, 0 ];
material.parameters{1}.sphere_radius        = sphere_radius(1);

%% Shape description for Fe
material.parameters{2}.name                 = 'Be';
material.parameters{2}.color_map            = [154 253 33]/255;

material.parameters{2}.sphere_centers       = [             0.25,            0.25,            0.25;
                                                     material.x3,     material.x3,     material.x3;
                                                     material.x3,     material.x3, 2-3*material.x3;
                                                     material.x3, 2-3*material.x3,     material.x3;
                                                 2-3*material.x3,     material.x3,     material.x3];
                                           
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 

%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};