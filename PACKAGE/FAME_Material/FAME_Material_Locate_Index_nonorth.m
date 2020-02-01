function [ org_idx ] = FAME_Material_Locate_Index_nonorth( Par_mesh, Par_lattice, Par_material )
% This programe return the indices corresponding to discrete points inside the material.
% Input: 
%    'grid_nums': must be a array which contains only 3 positive integers,
%           which stand for the discrete size for x-, y- and z- axis.
%    'lattice_type': must be a string array. It can be choosen as
%           'Simple_cubic', 'Face_centered_cubic' or 'Body_centered_cubic'
%    'lattice_constant': must be a positive real number.
%    'material_handle': it must be a function handle with the following
%           input and ouput forms:
%                index_in_material = material_handle(x,y,z,a1,a2,a3);
%           x,y,z are consider as the component of discrete points over the
%           parallel hexahedron spaned by a1, a2 and a3. For an example,
%           we may using following function handle to type spheres at
%           each corner of the parallel hexahedron corresponding to given a1, a2, a3.
%                    r = 0.4;
%                    sphere_handle = @(x,y,z,a1,a2,a3) find( (x).^2 + (y).^2 + (z).^2 < r^2 );
%           'material_handle' could be set as cell type parameter for two or more material. For example:
%                   gyroid_handle{1} = @(x,y,z,a1,a2,a3) FAME_Material_Locate_Handle_Gyroid(  x,  y,  z, a1, a2, a3, 1.1, Popt.domain.lattice_constant );
%                   gyroid_handle{2} = @(x,y,z,a1,a2,a3) FAME_Material_Locate_Handle_Gyroid( -x, -y, -z, a1, a2, a3, 1.1, Popt.domain.lattice_constant );
%    'display_material': defalt seting is 'off'. This function display the
%           material grid in a figure(12) if 'display_material' is put as a string
%           'on'.
% Output: the inside-material indices for Yee's grid, corresponding to these input.
    [ Origin_point_set] = FAME_Matrix_Grid_nonorth( Par_mesh.grid_num, Par_lattice.lattice_vec_a );                              

    if isfield(Par_material,'color_map') == 0
        Par_material.color_map = { [15 131 225]/255, [200,40,45]/255, [238,126,2]/255, [75,165,102]/255, rand(1,3), rand(1,3), rand(1,3), rand(1,3), rand(1,3), rand(1,3)};
    end
    
    %% Constructing the discretized grids
    [  org_idx, reshaped_point_set  ] = Material_Locate( Origin_point_set, Par_lattice.lattice_vec_a, Par_lattice.Omega, Par_lattice.lattice_constant.Permutation, Par_material.material_handle, Par_material.periodic_judge );

    %% Plot Original material grids
    if iscell(org_idx) == 1
        material_num = length(org_idx);   
    else
        material_num   = 1;
        org_idx        = {org_idx};
    end
    
    if strcmp( Par_material.display_material, 'on')
        figure(12); hold on
        for i = 1:material_num
            plot_material(reshaped_point_set, org_idx{i}, Par_lattice.Omega'*Par_lattice.lattice_vec_a, Par_lattice.lattice_constant, Par_material.color_map{i}, 'original');
        end
        title('Material grids in original domain')
    end     
end

function [ point_set_idx, reshaped_point_set ] = Material_Locate( point_set, lattice_vec_a, Omega, Permutation, material_handle, periodic_judge )
    point_set_orig       = point_set*Omega;
    lattice_vec_a_orig_P = Omega'*lattice_vec_a;
    I = eye(3); P = I(:,Permutation); invP = P';
    invPermutation = [find(invP(:,1)==1), find(invP(:,2)==1), find(invP(:,3)==1)];
    lattice_vec_a_orig   = lattice_vec_a_orig_P(:,invPermutation);

    invAP  = inv(lattice_vec_a_orig_P);
    coef   = point_set_orig*invAP';

    shift_1 = floor(coef(:,1));
    shift_2 = floor(coef(:,2));
    shift_3 = floor(coef(:,3));

    point_set_orig = point_set_orig - [shift_1,shift_2,shift_3]*lattice_vec_a_orig_P'; 

    % Compute 27 shift values
    if strcmp(periodic_judge,'on')
        e = ones(3,1);  s = [-1;0;1];
        shift = [kron(kron(e,e),s), kron(kron(e,s),e), kron(kron(s,e),e)];
        
        point_set_idx = cell(0);
        for i = 1:size(shift,1)
            shift_value    = shift(i,:)*lattice_vec_a_orig';
            point_set_idx_tmp  = material_handle(    point_set_orig(:,1)+shift_value(1),     point_set_orig(:,2)+shift_value(2),     point_set_orig(:,3)+shift_value(3), ...
                                                                lattice_vec_a_orig(:,1),                lattice_vec_a_orig(:,2),                lattice_vec_a_orig(:,3) );
            for j = 1:length(point_set_idx_tmp)
                if i == 1
                    point_set_idx{j} = [];
                end
                point_set_idx{j} = union( point_set_idx{j}, point_set_idx_tmp{j} );
            end
        end
    elseif strcmp(periodic_judge,'off')
        point_set_idx  = material_handle(     point_set_orig(:,1),     point_set_orig(:,2),     point_set_orig(:,3), ...
                                          lattice_vec_a_orig(:,1), lattice_vec_a_orig(:,2), lattice_vec_a_orig(:,3) );
    end

    reshaped_point_set = point_set_orig;
end


function plotmesh3D(mesh_x,mesh_y,mesh_z) 
    hold on
    for i = 1:length(mesh_z)
        [tmp_grid_x,tmp_grid_y] = meshgrid(mesh_x,mesh_y);
        mesh( tmp_grid_x, tmp_grid_y, mesh_z(i)*ones(size(tmp_grid_x)),'LineWidth',1 );
    end
    for i = 1:length(mesh_x)
        [tmp_grid_y,tmp_grid_z] = meshgrid(mesh_y,[mesh_z(1),mesh_z(end)] );
        mesh( mesh_x(i)*ones(size(tmp_grid_y)), tmp_grid_y, tmp_grid_z,'LineWidth',1  );
    end
    hidden off 
    colormap([0,0,0])
    axis equal
%     axis off
end

function plot_material(point_set, point_set_idx, lattice_vec_a, lattice_constant, color_map,lattice_vec_mode)
    plot3(point_set(point_set_idx,1), point_set(point_set_idx,2), point_set(point_set_idx,3), 'o',...
          'MarkerSize',2,...
          'MarkeredgeColor',color_map,...
          'MarkerfaceColor',color_map);
    if isempty(lattice_vec_a) == 0
        FAME_Plot_Parallelepiped([0,0,0],lattice_vec_a(:,1),lattice_vec_a(:,2),lattice_vec_a(:,3), lattice_constant,'color',lattice_vec_mode);
    end
    axis equal
    axis tight
%     axis off

    view([-15 0])
end