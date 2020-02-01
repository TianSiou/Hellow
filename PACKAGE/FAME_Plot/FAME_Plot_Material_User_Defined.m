function h = FAME_Plot_Material_User_Defined(grid_nums,a1,a2,a3,lattice_type,lattice_constant,Omega,Par_material,hax,plot_mode)
%     grid_nums = varargin{1}; 
%     a1 = varargin{2};  a2 = varargin{3};  a3 = varargin{4};
%     lattice_type = varargin{5};  lattice_constant = varargin{6};  Omega = varargin{7};  
%     switch nargin
%         case 10
%             function_handle_type = 'Sphere_and_Cylinder';
%             Par_material = varargin{8};  
%             hax = varargin{9};  
%             plot_mode = varargin{10};
%         case 11
%             function_handle_type = 'Isofunction';
%             isofunction = varargin{8};  
%             isovalue = varargin{9};  
%             hax = varargin{10};  
%             plot_mode = varargin{11};
%         otherwise
%             error('To many or less input arguments!');
%     end
    
    if isfield(Par_material, 'isofunction')
        function_handle_type = 'Isofunction';
        isofunction = Par_material.isofunction;
        isovalue    = Par_material.isovalue;
    else
        function_handle_type = 'Sphere_and_Cylinder';
        periodic_judge = Par_material.periodic_judge;
        parameters     = Par_material.data.parameters;
    end
    
    axes(hax); cla
    hold on
    %% Compute orginal lattice vectors
    lattice_vec_a        = [a1,a2,a3];
    lattice_vec_a_orig_P = Omega'*lattice_vec_a;
    I = eye(3); P = I(:,lattice_constant.Permutation); invP = P';
    invPermutation = [find(invP(:,1)==1), find(invP(:,2)==1), find(invP(:,3)==1)];
    lattice_vec_a_orig   = lattice_vec_a_orig_P(:,invPermutation);
    
    switch plot_mode
        case 'primitive_cell'
            invAP  = inv(lattice_vec_a_orig_P);
            Plot_Parallelepiped([0,0,0],lattice_vec_a_orig(:,1),lattice_vec_a_orig(:,2),lattice_vec_a_orig(:,3),'color');
            rate = 1.1;
            text( rate*lattice_vec_a_orig(1,1), rate*lattice_vec_a_orig(2,1), rate*lattice_vec_a_orig(3,1), '\bf{a}_1', 'color', 'black', 'FontSize', 14);
            text( rate*lattice_vec_a_orig(1,2), rate*lattice_vec_a_orig(2,2), rate*lattice_vec_a_orig(3,2), '\bf{a}_2', 'color', 'black', 'FontSize', 14);
            text( rate*lattice_vec_a_orig(1,3), rate*lattice_vec_a_orig(2,3), rate*lattice_vec_a_orig(3,3), '\bf{a}_3', 'color', 'black', 'FontSize', 14);
            title('Primitive cell')
        case 'unit_cell'
            convensional_lattice_vec_a = Conventional_Lattice_Vector( lattice_type, lattice_vec_a_orig );
            invAP  = inv(convensional_lattice_vec_a);
            Plot_Parallelepiped([0,0,0],convensional_lattice_vec_a(:,1),convensional_lattice_vec_a(:,2),convensional_lattice_vec_a(:,3),'no_color');
            Plot_Parallelepiped([0,0,0],lattice_vec_a_orig(:,1),lattice_vec_a_orig(:,2),lattice_vec_a_orig(:,3),'color');
            title('Unit cell')
        case 'computational_cell'
            computational_lattice_vec_a = Omega'*[lattice_vec_a(1,1), 0 0; 0, lattice_vec_a(2,2), 0; 0, 0, lattice_vec_a(3,3)];
            invAP  = inv(computational_lattice_vec_a);
            Plot_Parallelepiped([0,0,0],computational_lattice_vec_a(:,1),computational_lattice_vec_a(:,2),computational_lattice_vec_a(:,3),'no_color');
            Plot_Parallelepiped([0,0,0],lattice_vec_a_orig_P(:,1),lattice_vec_a_orig_P(:,2),lattice_vec_a_orig_P(:,3),'color');
            title('Computational cell')
        otherwise
            error('Invalid plot mode!');
    end

    %% Construct mesh grid
    if strcmp(plot_mode, 'primitive_cell')
        temp = [[0;0;0], lattice_vec_a_orig(:,1),lattice_vec_a_orig(:,2),lattice_vec_a_orig(:,3),lattice_vec_a_orig(:,1)+lattice_vec_a_orig(:,2),lattice_vec_a_orig(:,1)+lattice_vec_a_orig(:,3),lattice_vec_a_orig(:,2)+lattice_vec_a_orig(:,3), lattice_vec_a_orig(:,1)+lattice_vec_a_orig(:,2)+lattice_vec_a_orig(:,3)];
    elseif strcmp(plot_mode, 'unit_cell')
        temp = [[0;0;0], convensional_lattice_vec_a(:,1),convensional_lattice_vec_a(:,2),convensional_lattice_vec_a(:,3),convensional_lattice_vec_a(:,1)+convensional_lattice_vec_a(:,2),convensional_lattice_vec_a(:,1)+convensional_lattice_vec_a(:,3),convensional_lattice_vec_a(:,2)+convensional_lattice_vec_a(:,3), convensional_lattice_vec_a(:,1)+convensional_lattice_vec_a(:,2)+convensional_lattice_vec_a(:,3)];
    elseif strcmp(plot_mode, 'computational_cell')
        temp = [[0;0;0], computational_lattice_vec_a(:,1),computational_lattice_vec_a(:,2),computational_lattice_vec_a(:,3),computational_lattice_vec_a(:,1)+computational_lattice_vec_a(:,2),computational_lattice_vec_a(:,1)+computational_lattice_vec_a(:,3),computational_lattice_vec_a(:,2)+computational_lattice_vec_a(:,3), computational_lattice_vec_a(:,1)+computational_lattice_vec_a(:,2)+computational_lattice_vec_a(:,3)];
    end
    e    = 2./grid_nums;
    x_start = min(temp(1,:))-e(1);  x_end   = max(temp(1,:))+e(1);
    y_start = min(temp(2,:))-e(2);  y_end   = max(temp(2,:))+e(2);
    z_start = min(temp(3,:))-e(3);  z_end   = max(temp(3,:))+e(3);
    X = linspace(x_start,x_end,grid_nums(1));  Y = linspace(y_start,y_end,grid_nums(2));  Z = linspace(z_start,z_end,grid_nums(3));
    [X,Y,Z] = meshgrid(X,Y,Z);
    %% Compute inner cell indices
    coef            = [X(:),Y(:),Z(:)]*invAP';
    
    idx_outofcell = union( find(coef(:,1)<0 | coef(:,1)>1),...
                           union( find(coef(:,2)<0 | coef(:,2)>1),find(coef(:,3)<0 | coef(:,3)>1) ) );

    color_map_set = { [15 131 225]/255, [200,40,45]/255, [238,126,2]/255, [75,165,102]/255, rand(1,3), rand(1,3), rand(1,3), rand(1,3), rand(1,3), rand(1,3)};
    %% Plot isosurface    
    switch function_handle_type
        case 'Sphere_and_Cylinder'
            p_sphere   = cell(length(parameters),1);
            p_cylinder = cell(length(parameters),1);
            for i = 1:length(parameters)
                if isfield(parameters{i},'color_map')
                    color_map = parameters{i}.color_map;
                else
                    color_map = color_map_set{i};
                end
                sphere_centers = parameters{i}.sphere_centers;
                sphere_radius  = parameters{i}.sphere_radius;
                cylinder_bot_centers = parameters{i}.cylinder_bot_centers;
                cylinder_top_centers = parameters{i}.cylinder_top_centers;
                cylinder_radius      = parameters{i}.cylinder_radius;

                G_sphere   = Compute_MinDist_Sphere([X(:),Y(:),Z(:)], lattice_vec_a_orig, sphere_centers,periodic_judge);
                G_cylinder = Compute_MinDist_Cylinder([X(:),Y(:),Z(:)], lattice_vec_a_orig, cylinder_bot_centers, cylinder_top_centers,periodic_judge);

                G_sphere   = 1./G_sphere;     G_sphere(idx_outofcell)   = 1./sphere_radius(1);
                G_cylinder = 1./G_cylinder;   G_cylinder(idx_outofcell) = 1./cylinder_radius(1);

                p_sphere{i} = patch( isosurface(X,Y,Z,reshape(G_sphere,grid_nums(1), grid_nums(2), grid_nums(3)),1./sphere_radius(i)) );
                p_sphere{i}.FaceColor = color_map;
                p_sphere{i}.EdgeColor = 'none';

                p_cylinder{i} = patch( isosurface(X,Y,Z,reshape(G_cylinder, grid_nums(1), grid_nums(2), grid_nums(3)),1./cylinder_radius(i)) );
                p_cylinder{i}.FaceColor = color_map;
                p_cylinder{i}.EdgeColor = 'none';
            end
        case 'Isofunction'
            G = isofunction(X(:),Y(:),Z(:),lattice_vec_a_orig(:,1),lattice_vec_a_orig(:,2),lattice_vec_a_orig(:,3));
            color_map = color_map_set;
            p = cell(length(G),1);
            for i = 1:length(G)
                G{i}(idx_outofcell) = 0;

                p{i} = patch( isosurface(X,Y,Z,reshape(G{i},grid_nums(2), grid_nums(1), grid_nums(3)),isovalue(i)) );
                p{i}.FaceColor = color_map{i};
                p{i}.EdgeColor = 'none';
            end
    end
    
    camlight right
    rotate3d on
    axis equal
    axis off
    view([-38 17])
   
end

function G = Compute_MinDist_Sphere(point_set_orig, lattice_vec_a_orig, sphere_centers, periodic_judge)
    G = inf(size(point_set_orig,1),1);
    % Compute 27 shift values
    if strcmp(periodic_judge,'on')
        e = ones(3,1);  s = [-1;0;1];
        shift = [kron(kron(e,e),s), kron(kron(e,s),e), kron(kron(s,e),e)];
        for i = 1:size(sphere_centers,1)
            exact_sphere_centers = sphere_centers(i,:)*lattice_vec_a_orig';
            for j = 1:size(shift,1)
                shift_value    = shift(j,:)*lattice_vec_a_orig';
%                 G_tmp = sqrt( (point_set_orig + shift_value - exact_sphere_centers).^2 * ones(3,1) );
                G_tmp = sqrt( (point_set_orig(:,1) + shift_value(1) - exact_sphere_centers(1)).^2 + ...
                              (point_set_orig(:,2) + shift_value(2) - exact_sphere_centers(2)).^2 + ...
                              (point_set_orig(:,3) + shift_value(3) - exact_sphere_centers(3)).^2        );
                G = min(G,G_tmp);
            end
        end
    elseif strcmp(periodic_judge,'off')
        for i = 1:size(sphere_centers,1)
            exact_sphere_centers = sphere_centers(i,:)*lattice_vec_a_orig';
%             G_tmp = sqrt( (point_set_orig - exact_sphere_centers).^2 * ones(3,1) );
            G_tmp = sqrt( (point_set_orig(:,1) - exact_sphere_centers(1)).^2 + ...
                          (point_set_orig(:,2) - exact_sphere_centers(2)).^2 + ...
                          (point_set_orig(:,3) - exact_sphere_centers(3)).^2        );
            G = min(G,G_tmp);
        end
    end
end

function G =  Compute_MinDist_Cylinder(point_set_orig, lattice_vec_a_orig, cylinder_bot_centers, cylinder_top_centers,periodic_judge)
    Point_num = size(point_set_orig,1);
    G = inf(Point_num,1);
    if strcmp(periodic_judge,'on')
        e = ones(3,1);  s = [-1;0;1];
        shift = [kron(kron(e,e),s), kron(kron(e,s),e), kron(kron(s,e),e)];
        for i = 1:size(cylinder_bot_centers,1)
            exact_cylinder_bot_centers = cylinder_bot_centers(i,:)*lattice_vec_a_orig';
            exact_cylinder_top_centers = cylinder_top_centers(i,:)*lattice_vec_a_orig';
            Direct_vector = exact_cylinder_top_centers - exact_cylinder_bot_centers;
            Direct_vector = Direct_vector / norm( Direct_vector );
            for j = 1:size(shift,1)
                shift_value   = shift(j,:)*lattice_vec_a_orig';
                Point_bot_set = point_set_orig + shift_value - exact_cylinder_bot_centers;
                Point_top_set = point_set_orig + shift_value - exact_cylinder_top_centers;
                Point_line_dist = sqrt( ( Point_bot_set(:,2)*Direct_vector(3) - Point_bot_set(:,3)*Direct_vector(2) ).^2 + ...
                                        ( Point_bot_set(:,3)*Direct_vector(1) - Point_bot_set(:,1)*Direct_vector(3) ).^2 + ... 
                                        ( Point_bot_set(:,1)*Direct_vector(2) - Point_bot_set(:,2)*Direct_vector(1) ).^2        );
                Point_line_bot_proj =    Point_bot_set(:,1)*Direct_vector(1) + Point_bot_set(:,2)*Direct_vector(2) + Point_bot_set(:,3)*Direct_vector(3);
                Point_line_top_proj = - (Point_top_set(:,1)*Direct_vector(1) + Point_top_set(:,2)*Direct_vector(2) + Point_top_set(:,3)*Direct_vector(3) );
                Point_line_proj     = min( Point_line_bot_proj ,  Point_line_top_proj );
                
                Point_line_dist(Point_line_proj <= 0) = 100000;
                G = min(G,Point_line_dist);
            end
        end
    elseif strcmp(periodic_judge,'off')
        for i = 1:size(cylinder_bot_centers,1)
            exact_cylinder_bot_centers = cylinder_bot_centers(i,:)*lattice_vec_a_orig';
            exact_cylinder_top_centers = cylinder_top_centers(i,:)*lattice_vec_a_orig';
            Direct_vector = exact_cylinder_top_centers - exact_cylinder_bot_centers;
            Direct_vector = Direct_vector / norm( Direct_vector );
            Point_bot_set = point_set_orig - exact_cylinder_bot_centers;
            Point_top_set = point_set_orig - exact_cylinder_top_centers;
            Point_line_dist   = sqrt( ( Point_bot_set(:,2)*Direct_vector(3) - Point_bot_set(:,3)*Direct_vector(2) ).^2 + ...
                                      ( Point_bot_set(:,3)*Direct_vector(1) - Point_bot_set(:,1)*Direct_vector(3) ).^2 + ... 
                                      ( Point_bot_set(:,1)*Direct_vector(2) - Point_bot_set(:,2)*Direct_vector(1) ).^2        );
            Point_line_bot_proj =    Point_bot_set(:,1)*Direct_vector(1) + Point_bot_set(:,2)*Direct_vector(2) + Point_bot_set(:,3)*Direct_vector(3);
            Point_line_top_proj = - (Point_top_set(:,1)*Direct_vector(1) + Point_top_set(:,2)*Direct_vector(2) + Point_top_set(:,3)*Direct_vector(3) );
            Point_line_proj     = min( [ Point_line_bot_proj , Point_line_top_proj ],[],2 );

            Point_line_dist(Point_line_proj <= 0) = 100000;
            G = min(G,Point_line_dist);
        end
    end
end

function convensional_lattice_vec_a = Conventional_Lattice_Vector( lattice_type, lattice_vec_a_orig )
    switch lattice_type
        case {'simple_cubic','triclinic','primitive_monoclinic','primitive_tetragonal','primitive_othorhombic','hexagonal','rhombohedral'}
            convensional_lattice_vec_a = lattice_vec_a_orig;
        case {'face_centered_cubic','face_centered_orthorhombic'}
            convensional_lattice_vec_a(:,1) = -lattice_vec_a_orig(:,1) + lattice_vec_a_orig(:,2) + lattice_vec_a_orig(:,3);
            convensional_lattice_vec_a(:,2) =  lattice_vec_a_orig(:,1) - lattice_vec_a_orig(:,2) + lattice_vec_a_orig(:,3);
            convensional_lattice_vec_a(:,3) =  lattice_vec_a_orig(:,1) + lattice_vec_a_orig(:,2) - lattice_vec_a_orig(:,3);
        case {'body_centered_cubic','body_centered_tetragonal','body_centered_orthorhombic'}
            convensional_lattice_vec_a(:,1) = lattice_vec_a_orig(:,2) + lattice_vec_a_orig(:,3);
            convensional_lattice_vec_a(:,2) = lattice_vec_a_orig(:,1) + lattice_vec_a_orig(:,3);
            convensional_lattice_vec_a(:,3) = lattice_vec_a_orig(:,1) + lattice_vec_a_orig(:,2);
        case 'base_centered_monoclinic'
            convensional_lattice_vec_a(:,1) = lattice_vec_a_orig(:,1) - lattice_vec_a_orig(:,2);
            convensional_lattice_vec_a(:,2) = lattice_vec_a_orig(:,1) + lattice_vec_a_orig(:,2);
            convensional_lattice_vec_a(:,3) = lattice_vec_a_orig(:,3);
        case 'c_base_centered_orthorhombic'
            convensional_lattice_vec_a(:,1) = lattice_vec_a_orig(:,1) + lattice_vec_a_orig(:,2);
            convensional_lattice_vec_a(:,2) = lattice_vec_a_orig(:,2) - lattice_vec_a_orig(:,1);
            convensional_lattice_vec_a(:,3) = lattice_vec_a_orig(:,3);
        otherwise
            error('Invalid lattice type!');
    end
end

function Plot_Parallelepiped(start_pt,vec1,vec2,vec3,flag_color)
    vec1     = reshape(vec1,3,1);
    vec2     = reshape(vec2,3,1);
    vec3     = reshape(vec3,3,1);
    start_pt = reshape(start_pt,3,1);
    
    vertex = [  start_pt(1) start_pt(2) start_pt(3);
                [start_pt(1) start_pt(2) start_pt(3)] + vec1';
                [start_pt(1) start_pt(2) start_pt(3)] + vec1'+vec2';
                [start_pt(1) start_pt(2) start_pt(3)] + vec2';
                [start_pt(1) start_pt(2) start_pt(3)] + vec3';
                [start_pt(1) start_pt(2) start_pt(3)] + vec3'+vec1';
                [start_pt(1) start_pt(2) start_pt(3)] + vec3'+vec1'+vec2';
                [start_pt(1) start_pt(2) start_pt(3)] + vec3'+vec2';];
    edge = [1 2;
            2 3;
            3 4;
            4 1;
            1 5;
            2 6;
            3 7;
            4 8;
            5 6;
            6 7;
            7 8;
            8 5];
    for i = 1:size(edge,1)        
        idx1 = edge(i,1);
        idx2 = edge(i,2);
        path = [vertex(idx1,:);vertex(idx2,:)];
        plot3( path(:,1), path(:,2), path(:,3), 'k-')
        hold on
    end
    if strcmp(flag_color,'no_color') == 1
        plot3( [start_pt(1);start_pt(1)+vec1(1)'], [start_pt(2);start_pt(2)+vec1(2)'], [start_pt(3);start_pt(3)+vec1(3)'], 'k-',...
               [start_pt(1);start_pt(1)+vec2(1)'], [start_pt(2);start_pt(2)+vec2(2)'], [start_pt(3);start_pt(3)+vec2(3)'], 'k-',...
               [start_pt(1);start_pt(1)+vec3(1)'], [start_pt(2);start_pt(2)+vec3(2)'], [start_pt(3);start_pt(3)+vec3(3)'], 'k-');
    elseif strcmp(flag_color,'color') == 1
        mArrow3([start_pt(1) start_pt(2) start_pt(3)],[start_pt(1)+vec1(1)' start_pt(2)+vec1(2)' start_pt(3)+vec1(3)'], 'facealpha', 1, 'color', 'red', 'stemWidth', 0.02, 'tipWidth', 0.04); 
        mArrow3([start_pt(1) start_pt(2) start_pt(3)],[start_pt(1)+vec2(1)' start_pt(2)+vec2(2)' start_pt(3)+vec2(3)'], 'facealpha', 1, 'color', 'green', 'stemWidth', 0.02, 'tipWidth', 0.04); 
        mArrow3([start_pt(1) start_pt(2) start_pt(3)],[start_pt(1)+vec3(1)' start_pt(2)+vec3(2)' start_pt(3)+vec3(3)'], 'facealpha', 1, 'color', 'blue', 'stemWidth', 0.02, 'tipWidth', 0.04); 
%         plot3( [start_pt(1);start_pt(1)+vec1(1)'], [start_pt(2);start_pt(2)+vec1(2)'], [start_pt(3);start_pt(3)+vec1(3)'], 'r-',...
%                [start_pt(1);start_pt(1)+vec2(1)'], [start_pt(2);start_pt(2)+vec2(2)'], [start_pt(3);start_pt(3)+vec2(3)'], 'g-',...
%                [start_pt(1);start_pt(1)+vec3(1)'], [start_pt(2);start_pt(2)+vec3(2)'], [start_pt(3);start_pt(3)+vec3(3)'], 'b-');
    end
end