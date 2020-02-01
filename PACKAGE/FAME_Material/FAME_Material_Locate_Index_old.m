function [ ele_x_idx, ele_y_idx, ele_z_idx, mag_x_idx, mag_y_idx, mag_z_idx, org_idx ] = FAME_Material_Locate_Index_old( Par_mesh, Par_lattice, Par_material )
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
    [ Electric_x_point_set, Electric_y_point_set, Electric_z_point_set, Magnetic_x_point_set, Magnetic_y_point_set, Magnetic_z_point_set,Origin_point_set] = ...
        FAME_Matrix_Grid( Par_mesh.grid_num, Par_mesh.mesh_len );                              


    if isfield(Par_material,'color_map') == 0
        Par_material.color_map = { [15 131 225]/255, [200,40,45]/255, [238,126,2]/255, [75,165,102]/255, rand(1,3), rand(1,3), rand(1,3), rand(1,3), rand(1,3), rand(1,3)};
    end
    
    %% Constructing the discretized grids on edge and face
    % Determining whether all the points in the material or not
%     if strcmp(displar_material,'on')
%         fig_grid = figure(10);clf
%         set(fig_grid,'name','Material grids');
%         hax_grid = axes(fig_grid);
%     end
    
    [ ele_x_idx ] = FAME_Material_Locate( Electric_x_point_set, Par_lattice.lattice_vec_a, Par_lattice.lattice_constant, Par_lattice.Omega, Par_material.material_handle, 'off' );
    [ ele_y_idx ] = FAME_Material_Locate( Electric_y_point_set, Par_lattice.lattice_vec_a, Par_lattice.lattice_constant, Par_lattice.Omega, Par_material.material_handle, 'off' );
    [ ele_z_idx ] = FAME_Material_Locate( Electric_z_point_set, Par_lattice.lattice_vec_a, Par_lattice.lattice_constant, Par_lattice.Omega, Par_material.material_handle, 'off' );

    [ mag_x_idx ] = FAME_Material_Locate( Magnetic_x_point_set, Par_lattice.lattice_vec_a, Par_lattice.lattice_constant, Par_lattice.Omega, Par_material.material_handle, 'off' );
    [ mag_y_idx ] = FAME_Material_Locate( Magnetic_y_point_set, Par_lattice.lattice_vec_a, Par_lattice.lattice_constant, Par_lattice.Omega, Par_material.material_handle, 'off' );
    [ mag_z_idx ] = FAME_Material_Locate( Magnetic_z_point_set, Par_lattice.lattice_vec_a, Par_lattice.lattice_constant, Par_lattice.Omega, Par_material.material_handle, 'off' );

    [  org_idx, reshaped_point_set  ] = FAME_Material_Locate( Origin_point_set, Par_lattice.lattice_vec_a, Par_lattice.lattice_constant, Par_lattice.Omega, Par_material.material_handle, 'off' );

    %% Plot Original material grids
    if iscell(org_idx) == 1
        material_num = length(org_idx);   
    else
        material_num   = 1;
        org_idx        = {org_idx};
        ele_x_idx      = {ele_x_idx};
        ele_y_idx      = {ele_y_idx};
        ele_z_idx      = {ele_z_idx};
        mag_x_idx      = {mag_x_idx};
        mag_y_idx      = {mag_y_idx};
        mag_z_idx      = {mag_z_idx};
    end

    if strcmp( Par_material.display_material, 'on')
        figure(12); subplot(1,2,1); hold on
        for i = 1:material_num
            plot_material(reshaped_point_set, org_idx{i}, Par_lattice.Omega'*Par_lattice.lattice_vec_a, Par_lattice.lattice_constant, Par_material.color_map{i}, 'original');
%             plot_material(reshaped_point_set, org_idx{i}, Omega\lattice_vec_a, lattice_constant, color_map{i}, 'original');
        end
        title('Material grids in original domain')
    end
% X = linspace(0,1,60);
% Y = linspace(0,1/sqrt(2),60);
% Z = linspace(0,1/sqrt(2),60);
% [X,Y,Z] = meshgrid(X,Y,Z);
% Origin_point_set = [X(:),Y(:),Z(:)];
    %% Plot Computational material grids 
    if strcmp( Par_material.display_material, 'on')
        figure(12); subplot(1,2,2); hold on
        plotmesh3D([0,Par_lattice.lattice_vec_a(1,1)],...
                   [0,Par_lattice.lattice_vec_a(2,2)],...
                   [0,Par_lattice.lattice_vec_a(3,3)]);
        for i = 1:material_num
                plot_material(Origin_point_set, org_idx{i}, Par_lattice.lattice_vec_a, Par_lattice.lattice_constant, Par_material.color_map{i}, 'computational');
        end
        title('Material grids in computational domain')
    end        
end

function [ point_set_idx, reshaped_point_set ] = FAME_Material_Locate( point_set, lattice_vec_a, lattice_constant, Omega, material_handle, display_mesh )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_xz = cross( lattice_vec_a(:,3),lattice_vec_a(:,2) );
temp_idx_xz_l = find( test_xz(3)*point_set(:,3) > -test_xz(1)*point_set(:,1) - test_xz(2)*point_set(:,2) + 1e-8  );
% size(temp_idx_xz_l)
coef_xz_l     = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_idx_xz_r = find( test_xz(3)*point_set(:,3) < -test_xz(1)*( point_set(:,1)-lattice_vec_a(1,1) ) - test_xz(2)*(point_set(:,2)-lattice_vec_a(2,1)) + 1e-8 );
% size(temp_idx_xz_r)
coef_xz_r     = -1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_yz = cross( lattice_vec_a(:,1),lattice_vec_a(:,3) );
temp_idx_yz_l = find( test_yz(3)*point_set(:,3) > -test_yz(1)*point_set(:,1) - test_yz(2)*point_set(:,2) + 1e-8  );
% size(temp_idx_yz_l)
coef_yz_l     = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_idx_yz_r = find( test_yz(3)*point_set(:,3) < -test_yz(1)*( point_set(:,1)-lattice_vec_a(1,2) ) - test_yz(2)*(point_set(:,2)-lattice_vec_a(2,2)) + 1e-8  );
% size(temp_idx_yz_r)
coef_yz_r     = -1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx_l_1    =   setdiff(setdiff( temp_idx_xz_l, temp_idx_yz_l ),temp_idx_yz_r);
idx_r_1    =   setdiff(setdiff( temp_idx_xz_r, temp_idx_yz_l ),temp_idx_yz_r);
idx_l_2    =   setdiff(setdiff( temp_idx_yz_l, temp_idx_xz_l ),temp_idx_xz_r);
idx_r_2    =   setdiff(setdiff( temp_idx_yz_r, temp_idx_xz_l ),temp_idx_xz_r);
idx_ll_3   =   intersect( temp_idx_xz_l, temp_idx_yz_l );
idx_lr_3   =   intersect( temp_idx_xz_l, temp_idx_yz_r );
idx_rl_3   =   intersect( temp_idx_xz_r, temp_idx_yz_l );
idx_rr_3   =   intersect( temp_idx_xz_r, temp_idx_yz_r );
% temp_res   =   setdiff( 1:size(point_set,1), union( idx_l_1, union(idx_r_1,union(idx_l_2,union(idx_r_2,union(idx_ll_3,union(idx_lr_3,union(idx_rl_3,idx_rr_3))))))));
idx_1      =   union(idx_l_1,idx_r_1) ;
idx_2      =   union(idx_l_2,idx_r_2) ;
idx_3      =   union(idx_ll_3,union(idx_lr_3,union(idx_rl_3,idx_rr_3 ))) ;
temp_res   =   setdiff( 1:size(point_set,1), union( idx_1, union(idx_2,idx_3)) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
point_set_temp     = point_set;
% only shift along a1
point_set(idx_l_1,1) = point_set(idx_l_1,1) + coef_xz_l*lattice_vec_a(1,1);
point_set(idx_l_1,2) = point_set(idx_l_1,2) + coef_xz_l*lattice_vec_a(2,1);
point_set(idx_r_1,1) = point_set(idx_r_1,1) + coef_xz_r*lattice_vec_a(1,1);
point_set(idx_r_1,2) = point_set(idx_r_1,2) + coef_xz_r*lattice_vec_a(2,1);
% only shift along a2
point_set(idx_l_2,1) = point_set(idx_l_2,1) + coef_yz_l*lattice_vec_a(1,2);
point_set(idx_l_2,2) = point_set(idx_l_2,2) + coef_yz_l*lattice_vec_a(2,2);
point_set(idx_r_2,1) = point_set(idx_r_2,1) + coef_yz_r*lattice_vec_a(1,2);
point_set(idx_r_2,2) = point_set(idx_r_2,2) + coef_yz_r*lattice_vec_a(2,2);
% only along a1+a2
point_set(idx_ll_3,1) = point_set(idx_ll_3,1) + coef_xz_l*lattice_vec_a(1,1) + coef_yz_l*lattice_vec_a(1,2);
point_set(idx_ll_3,2) = point_set(idx_ll_3,2) + coef_xz_l*lattice_vec_a(2,1) + coef_yz_l*lattice_vec_a(2,2);
point_set(idx_lr_3,1) = point_set(idx_lr_3,1) + coef_xz_l*lattice_vec_a(1,1) + coef_yz_r*lattice_vec_a(1,2);
point_set(idx_lr_3,2) = point_set(idx_lr_3,2) + coef_xz_l*lattice_vec_a(2,1) + coef_yz_r*lattice_vec_a(2,2);
point_set(idx_rl_3,1) = point_set(idx_rl_3,1) + coef_xz_r*lattice_vec_a(1,1) + coef_yz_l*lattice_vec_a(1,2);
point_set(idx_rl_3,2) = point_set(idx_rl_3,2) + coef_xz_r*lattice_vec_a(2,1) + coef_yz_l*lattice_vec_a(2,2);
point_set(idx_rr_3,1) = point_set(idx_rr_3,1) + coef_xz_r*lattice_vec_a(1,1) + coef_yz_r*lattice_vec_a(1,2);
point_set(idx_rr_3,2) = point_set(idx_rr_3,2) + coef_xz_r*lattice_vec_a(2,1) + coef_yz_r*lattice_vec_a(2,2);
% point_set(idx_3,2) = point_set(idx_3,2) + coef_yz*lattice_vec_a(2,2);

point_set_orig     = point_set*Omega;
% point_set_orig     = (Omega\point_set')';
lattice_vec_a_orig = Omega'*lattice_vec_a;
% P = eye(3); P = P(:,lattice_constant.Permutation);
% lattice_vec_a_orig = Omega\lattice_vec_a*P'

%     figure(14); hold on
%     plot3([0,lattice_vec_a_orig(1,1)],[0,lattice_vec_a_orig(2,1)],[0,lattice_vec_a_orig(3,1)],'r-')
%     plot3([0,lattice_vec_a_orig(1,2)],[0,lattice_vec_a_orig(2,2)],[0,lattice_vec_a_orig(3,2)],'g-')
%     plot3([0,lattice_vec_a_orig(1,3)],[0,lattice_vec_a_orig(2,3)],[0,lattice_vec_a_orig(3,3)],'b-')
%     plot3(point_set_orig(:,1),point_set_orig(:,2),point_set_orig(:,3),'k.');
%     axis equal
    
% X = linspace(0,1,60);
% Y = linspace(0,1/sqrt(2),60);
% Z = linspace(0,1/sqrt(2),60);
% [X,Y,Z] = meshgrid(X,Y,Z);
% point_set_orig = [X(:),Y(:),Z(:)];

point_set_idx      = material_handle(    point_set_orig(:,1),     point_set_orig(:,2),     point_set_orig(:,3), ...
                                     lattice_vec_a_orig(:,1), lattice_vec_a_orig(:,2), lattice_vec_a_orig(:,3) );
%     figure(15); hold on
%     plot3([0,lattice_vec_a_orig(1,1)],[0,lattice_vec_a_orig(2,1)],[0,lattice_vec_a_orig(3,1)],'r-')
%     plot3([0,lattice_vec_a_orig(1,2)],[0,lattice_vec_a_orig(2,2)],[0,lattice_vec_a_orig(3,2)],'g-')
%     plot3([0,lattice_vec_a_orig(1,3)],[0,lattice_vec_a_orig(2,3)],[0,lattice_vec_a_orig(3,3)],'b-')
%     plot3(point_set_orig(point_set_idx,1),point_set_orig(point_set_idx,2),point_set_orig(point_set_idx,3),'k.');
%     axis equal
                                 
                                 
reshaped_point_set = point_set_orig;                                 
 
%% Plot computational mesh and reshaped mesh
if strcmp( display_mesh, 'on')
    figure(11)
    subplot(1,2,1); hold on
    FAME_Plot_Parallelepiped([0,0,0],lattice_vec_a(:,1),lattice_vec_a(:,2),lattice_vec_a(:,3),lattice_constant,'color', 'computational')
    plot3(point_set_temp(temp_res,1),point_set_temp(temp_res,2),point_set_temp(temp_res,3),'k.',...
          point_set_temp(idx_1,1),point_set_temp(idx_1,2),point_set_temp(idx_1,3),'r.',...
          point_set_temp(idx_2,1),point_set_temp(idx_2,2),point_set_temp(idx_2,3),'g.',...
          point_set_temp(idx_3,1),point_set_temp(idx_3,2),point_set_temp(idx_3,3),'b.')
    title('Computational mesh')
    axis equal
    axis tight
    
    subplot(1,2,2); hold on
    FAME_Plot_Parallelepiped([0,0,0],lattice_vec_a(:,1),lattice_vec_a(:,2),lattice_vec_a(:,3),lattice_constant,'color','computational')
    plot3(point_set(temp_res,1),point_set(temp_res,2),point_set(temp_res,3),'k.',...
          point_set(idx_1,1),point_set(idx_1,2),point_set(idx_1,3),'r.',...
          point_set(idx_2,1),point_set(idx_2,2),point_set(idx_2,3),'g.',...
          point_set(idx_3,1),point_set(idx_3,2),point_set(idx_3,3),'b.');
    title('Reshaped computational mesh')
    axis equal
    axis tight
end
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