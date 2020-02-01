function material = No000_Primitive_Tetragonal_testDemo_cube(sphere_radius, cylinder_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

material.Prototype    = 'No000_Primitive_Tetragonal_testDemo_cube';
material.material_num = 1;
ncell = 27;
N     = 3;

material.lattice_type = 'simple_cubic';
% material.lattice_type = 'primitive_tetragonal';

material.lattice_constant.a      = N;
% material.lattice_constant.b      = ncell;
% material.lattice_constant.c      =  ncell;

[sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof( sphere_radius, cylinder_radius, material.material_num );
%% Shape description for Po
material.parameters{1}.name                 = 'Po';
material.parameters{1}.color_map            = [168 109 57]/255;

for ii = 1 : N
   for jj = 1 : N
       for kk = 1 : N
            material.parameters{1}.sphere_centers( N*N*(ii-1)+N*(jj-1)+kk,:) = [1/(2*N)+ (kk-1)/N, 1/(2*N)+ (jj-1)/N, 1/(2*N)+ (ii-1)/N];
       end
   end
end
if N ~= 1
    material.parameters{1}.sphere_centers( (ncell+1)/2 , : ) = [];
end
material.parameters{1}.sphere_radius        = sphere_radius(1)*ones(1,size(material.parameters{1}.sphere_centers,1)); 

%% Shape description for Po (Cylinder)
for ii = 1 : N
    for jj = 1 : N
        for kk = 1 : N
material.parameters{1}.cylinder_bot_centers(1 + 3*(kk-1) + 3*N*N*(ii-1)+3*N*(jj-1): 3*kk+ 3*N*N*(ii-1)+3*N*(jj-1),:)       = [ 1/(2*N) + (kk-1)/N, 1/(2*N) + (jj-1)/N, (ii-1)/N ;
                                                                                                                        (kk-1)/N, 1/(2*N) + (jj -1)/N, 1/(2*N) + (ii -1)/N;
                                                                                                                         1/(2*N) + (kk-1)/N, (jj -1)/N ,1/(2*N)+(ii -1)/N];                                          
material.parameters{1}.cylinder_top_centers(1 + 3*(kk-1) + 3*N*N*(ii-1)+3*N*(jj-1): 3*kk+ 3*N*N*(ii-1)+3*N*(jj-1),:)        = [1/(2*N) + (kk-1)/N, 1/(2*N) + (jj-1)/N, ii/N ;
                                                                                                                        kk/N, 1/(2*N) + (jj -1)/N, 1/(2*N) + (ii -1)/N;
                                                                                                                         1/(2*N) + (kk-1)/N, jj/N ,1/(2*N)+(ii -1)/N];
        end
    end
end
if N ~= 1
material.parameters{1}.cylinder_bot_centers( (ncell+1)/2 *3 -2 : (ncell+1)/2 *3  , :) = [];
material.parameters{1}.cylinder_top_centers( (ncell+1)/2 *3 -2: (ncell+1)/2 *3  , :) = [];
end
material.parameters{1}.cylinder_radius           = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_top_centers,1)); 



% material.parameters{1}.cylinder_top_centers(N*N*(ii-1)+N*(jj-1)+kk,:)        = [1/(2*N) + (kk-1)/N, 1/(2*N) + (jj-1)/N, ii/N ];
% material.parameters{1}.cylinder_bot_centers(1 + 3*ii : 3*(ii+1), :)       = [ 0, 0, (ii ) / sum(n_supercell);
%                                                       0, 0, (ii ) / sum(n_supercell);
%                                                       0, 0, (ii ) / sum(n_supercell)];                                               
% material.parameters{1}.cylinder_top_centers(1 + 3*ii : 3*(ii+1), :)        = [ 1, 0, (ii) / sum(n_supercell);
%                                                       0, 1, (ii ) / sum(n_supercell);
%                                                       0, 0, ii / sum(n_supercell)];
% material.parameters{1}.cylinder_radius           = cylinder_radius(1)*ones(1,size(material.parameters{1}.cylinder_top_centers,1));
%% Collect color map
material.color_map = {material.parameters{1}.color_map};