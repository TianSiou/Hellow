function material = No000_Primitive_Tetragonal_testDemo_sup(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'No000_Primitive_Tetragonal_testDemo_sup';
material.material_num = 1;
n_supercell = [2, 2];

material.lattice_type = 'primitive_tetragonal';

material.lattice_constant.a      = 1;
material.lattice_constant.b      = 1;
material.lattice_constant.c      =  sum(n_supercell);

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Po
material.parameters{1}.name                 = 'Po';
material.parameters{1}.color_map            = [168 109 57]/255;
for ii = 1 : sum(n_supercell)
material.parameters{1}.sphere_centers(ii, :)       = [ 0, 0, (ii-1) / sum(n_supercell)];                              
end
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 

%% Shape description for Po (Cylinder)
for ii = 1 : sum(n_supercell) / 2
material.parameters{1}.cylinder_bot_centers(1 + 3*(ii-1) : 3*ii, :)       = [ 0, 0, (ii - 1) / sum(n_supercell);
                                                      0, 0, (ii - 1) / sum(n_supercell);
                                                      0, 0, (ii - 1) / sum(n_supercell)];                                               
material.parameters{1}.cylinder_top_centers(1 + 3*(ii-1) : 3*ii, :)        = [ 1, 0, (ii - 1) / sum(n_supercell);
                                                      0, 1, (ii - 1) / sum(n_supercell);
                                                      0, 0, ii / sum(n_supercell)];
material.parameters{1}.cylinder_radius           = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_top_centers,1)); 
end
material.parameters{1}.cylinder_bot_centers(1 + 3*ii : 3*(ii+1), :)       = [ 0, 0, (ii ) / sum(n_supercell);
                                                      0, 0, (ii ) / sum(n_supercell);
                                                      0, 0, (ii ) / sum(n_supercell)];                                               
material.parameters{1}.cylinder_top_centers(1 + 3*ii : 3*(ii+1), :)        = [ 1, 0, (ii) / sum(n_supercell);
                                                      0, 1, (ii ) / sum(n_supercell);
                                                      0, 0, ii / sum(n_supercell)];
material.parameters{1}.cylinder_radius           = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_top_centers,1));
%% Collect color map
material.color_map = {material.parameters{1}.color_map};