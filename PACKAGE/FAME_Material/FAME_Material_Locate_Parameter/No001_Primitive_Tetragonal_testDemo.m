function material = No001_Primitive_Tetragonal_testDemo(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'No001_Primitive_Tetragonal_testDemo';
material.material_num = 1;
ncell = [2, 2];
% if ncell == 1
% material.lattice_type = 'simple_cubic';
% else
material.lattice_type = 'primitive_tetragonal';
% end

material.lattice_constant.a      = 1;
material.lattice_constant.b      = 1;
material.lattice_constant.c      = sum(ncell);

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Po
material.parameters{1}.name                 = 'Po';
material.parameters{1}.color_map            = [168 109 57]/255;
for ii = 1 : sum(ncell)/2
material.parameters{1}.sphere_centers(ii, :)       = [ 1/2, 1/2, ( (ii - 1) + 1/2 )/sum(ncell)];                              
end
material.parameters{1}.sphere_radius       = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 
% material.parameters{1}.sphere_radius(1,   1:ii/2)        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 
% material.parameters{1}.sphere_radius(1,ii/2+1 : ii)        = 0.8 .*ones(1,ii/2); 
%% Shape description for Po (Cylinder)
for ii = 1 : sum(ncell)/2
material.parameters{1}.cylinder_bot_centers(1 + 3*(ii-1) : 3*ii, :)       = [ 0, 1/2, ( (ii - 1) + 1/2 )/ sum(ncell);
                                                      1/2, 0,  ( (ii - 1) + 1/2 )/ sum(ncell);
                                                      1/2, 1/2, (ii - 1) / sum(ncell);
                                                      ];                                               
material.parameters{1}.cylinder_top_centers(1 + 3*(ii-1) : 3*ii, :)        = [ 1, 1/2, ( (ii - 1) + 1/2 )/ sum(ncell)  ;
                                                      1/2, 1, ( (ii - 1) + 1/2 )/ sum(ncell);
                                                      1/2, 1/2, ii / sum(ncell)];
material.parameters{1}.cylinder_radius           = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_top_centers,1)); 
end
for i =1 : sum(ncell)/2 
material.parameters{1}.cylinder_bot_centers(3*sum(ncell)/2 + i , :)       =    [  1/2, 1/2, (ii - 1) / sum(ncell) *2];
material.parameters{1}.cylinder_top_centers(3*sum(ncell)/2 + i , :)       =    [1/2, 1/2, ii / sum(ncell) *2] ;
% material.parameters{1}.cylinder_bot_centers( (3*sum(ncell)/2) +ii, :)       = [  1/2, 1/2, (ii - 1) / sum(ncell)];                                            
% material.parameters{1}.cylinder_top_centers((3*sum(ncell)/2) +ii , :)        = [1/2, 1/2, ii / sum(ncell)] ;
% material.parameters{1}.cylinder_radius                        = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_top_centers,1));
% material.parameters{1,1}.cylinder_radius(:, size(material.parameters{1}.cylinder_top_centers,1)+i-1)           = 0.8 ;
end
material.parameters{1,1}.cylinder_radius(:, size(material.parameters{1}.cylinder_top_centers,1)-ii+1:size(material.parameters{1}.cylinder_top_centers,1))           = 0.8 .*ones(1, ii);
%% Collect color map
material.color_map = {material.parameters{1}.color_map};