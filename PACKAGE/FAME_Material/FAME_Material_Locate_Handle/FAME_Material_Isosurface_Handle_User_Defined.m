function fun_iso = FAME_Material_Isosurface_Handle_User_Defined( material )
    for i = 1:length(material.parameters)
        fun_iso{i} = @(x,y,z,a1,a2,a3) Dist_Combine(x,y,z,a1,a2,a3,...
            material.parameters{i}.sphere_centers, material.parameters{i}.sphere_radius,...
            material.parameters{i}.cylinder_bot_centers, material.parameters{i}.cylinder_top_centers, material.parameters{i}.cylinder_radius);
    end
end

function dist = Dist_Combine(x,y,z,a1,a2,a3,sphere_centers,sphere_radius,cylinder_bot_centers,cylinder_top_centers,cylinder_radius)
    x = x(:);
    y = y(:);
    z = z(:);
%     shift = [0,0,0;
%               0,  1,  0];
    shift = [ -1, -1, -1; 
               0, -1, -1;
               1, -1, -1;
              -1,  0, -1; 
               0,  0, -1;
               1,  0, -1;
              -1,  1, -1; 
               0,  1, -1;
               1,  1, -1; 
              -1, -1,  0; 
               0, -1,  0;
               1, -1,  0;
              -1,  0,  0; 
               0,  0,  0;
               1,  0,  0;
              -1,  1,  0; 
               0,  1,  0;
               1,  1,  0;
              -1, -1,  1; 
               0, -1,  1;
               1, -1,  1;
              -1,  0,  1; 
               0,  0,  1;
               1,  0,  1;
              -1,  1,  1; 
               0,  1,  1;
               1,  1,  1];

%     dist_shpere = zeros(size(x))';
    center_shift = ( sphere_centers(1,:) + shift(1,:) )*[a1,a2,a3]';
    dist_shpere = Dist_Sphere([x,y,z],sphere_radius(1), center_shift);
    dist_shpere = dist_shpere';
    for j = 1:size(sphere_centers,1)
        for k=2:size(shift,1)
            center_shift = ( sphere_centers(j,:) + shift(k,:) )*[a1,a2,a3]';
            temp = Dist_Sphere([x,y,z],sphere_radius(j), center_shift);
            dist_shpere = min(dist_shpere, temp');
        end
    end

    bot_center_shift = ( cylinder_bot_centers(1,:) + shift(1,:) )*[a1,a2,a3]';
    top_center_shift = ( cylinder_top_centers(1,:) + shift(1,:) )*[a1,a2,a3]';
    dist_cylinder    = Dist_Cylinder([x,y,z],cylinder_radius(1), bot_center_shift, top_center_shift);
    dist_cylinder    = dist_cylinder';
    for j = 1:size(cylinder_bot_centers,1)
        for k = 2:size(shift,1)
            bot_center_shift = ( cylinder_bot_centers(j,:) + shift(k,:) )*[a1,a2,a3]';
            top_center_shift = ( cylinder_top_centers(j,:) + shift(k,:) )*[a1,a2,a3]';
            temp = Dist_Cylinder([x,y,z],cylinder_radius(j), bot_center_shift, top_center_shift);
            dist_cylinder = min(dist_cylinder, temp');
        end
    end
    dist = [dist_shpere;dist_cylinder]';
end

function dist = Dist_Sphere(Point_set, sphere_radius, sphere_center)
    dist = sqrt( (Point_set(:,1)-sphere_center(1)).^2 + ...
                 (Point_set(:,2)-sphere_center(2)).^2 + ...
                 (Point_set(:,3)-sphere_center(3)).^2);
%     dist( dist > sphere_radius ) = 0;
end

function dist = Dist_Cylinder(Point_set, cylinder_radius, bot_center, top_center)
    % Counting the number of points
    Point_num = size(Point_set,1);
    % Shifting the points
    Point_bot_set = Point_set - kron( ones( Point_num , 1 ) , bot_center );
    Point_top_set = Point_set - kron( ones( Point_num , 1 ) , top_center );
    % Computing the direction of cylinder
    Direct_vector = top_center - bot_center;
    % Normalizing the direction vector
    Direct_vector =       Direct_vector / norm( Direct_vector );
    % Constructing the array for vectorization computing 
    Direct_vector_array = kron( ones( Point_num , 1 ) , Direct_vector );
    % Computing the distance for each point to the direction vector
    dist   = sqrt(cross( Point_bot_set , Direct_vector_array ).^2 * ones(3,1));
    % Computing norm of the projection vectors 
    Point_line_bot_proj =  Point_bot_set*Direct_vector.';
    Point_line_top_proj = -Point_top_set*Direct_vector.';
    Point_line_proj     = min( [ Point_line_bot_proj , Point_line_top_proj ],[],2 );
    % Determining whether Point_set in a cylinder or not
    Point_idx_proj = find( Point_line_proj >= 0 );
%     dist( dist > cylinder_radius^2 ) = 0;
%     dist( Point_idx_proj ) = inf;
end