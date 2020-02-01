function material = No215_SC_Fe4C(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No215_SC_Fe4C
%
% When material.x2 = 1/4, the iron atoms are at the positions of the face-centered cubic lattice.
% In Fe4C, material.x2 is about 0.265. (that is, the lattice constant a is about 4*0.265)
% More detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/AB4_cP5_215_a_e.html
%
% Edit at 2017/6/17 By Jia-Wei Lin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'Fe4C';
material.material_num = 2;

material.lattice_type = 'simple_cubic';
material.website      = 'http://www.aflowlib.org/CrystalDatabase/AB4_cP5_215_a_e.html';
material.x2           = 0.265;

material.lattice_constant.a = 3.878;
[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for C
material.parameters{1}.name                 = 'C';
material.parameters{1}.color_map            = [177 173 165]/255;
material.parameters{1}.sphere_centers       = [ 0, 0, 0 ];
material.parameters{1}.sphere_radius        = sphere_radius(1);

material.parameters{1}.cylinder_bot_centers = [        0,         0,         0;
                                                       1,         1,         0; 
                                                       1,         0,         1;
                                                       0,         1,         1];
material.parameters{1}.cylinder_top_centers = [   .5*material.x2 ,     .5*material.x2 ,     .5*material.x2 ; 
                                               .5*(2-material.x2),  .5*(2-material.x2),     .5*material.x2 ; 
                                               .5*(2-material.x2),     .5*material.x2 ,  .5*(2-material.x2);
                                                  .5*material.x2 ,  .5*(2-material.x2),  .5*(2-material.x2)];                                              
material.parameters{1}.cylinder_radius        = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_bot_centers,1)); 
%% Shape description for Fe
material.parameters{2}.name                 = 'Fe';
material.parameters{2}.color_map            = [241 130 82]/255;
material.parameters{2}.sphere_centers       = [   material.x2,   material.x2,   material.x2;
                                                1-material.x2, 1-material.x2,   material.x2;
                                                1-material.x2,   material.x2, 1-material.x2;
                                                  material.x2, 1-material.x2, 1-material.x2];                                            
material.parameters{2}.sphere_radius        = sphere_radius(2)*ones(1,size(material.parameters{2}.sphere_centers,1)); 
material.parameters{2}.cylinder_bot_centers = [    .5*material.x2 ,     .5*material.x2 ,     .5*material.x2 ; 
                                                .5*(2-material.x2),  .5*(2-material.x2),     .5*material.x2 ; 
                                                .5*(2-material.x2),     .5*material.x2 ,  .5*(2-material.x2);
                                                   .5*material.x2 ,  .5*(2-material.x2),  .5*(2-material.x2)];
material.parameters{2}.cylinder_top_centers = material.parameters{2}.sphere_centers;                                           
material.parameters{2}.cylinder_radius        = cylinder_radius(2)*ones(1,size(material.parameters{2}.cylinder_bot_centers,1)); 
%% Collect color map
material.color_map = {material.parameters{1}.color_map, material.parameters{2}.color_map};