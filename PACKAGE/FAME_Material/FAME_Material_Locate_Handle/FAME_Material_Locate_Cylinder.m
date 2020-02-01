function [ Point_idx ] = FAME_Material_Locate_Cylinder( Point_set, Cylinder_radius, Cylinder_bot_center, Cylinder_top_center )
%% This function is used to determine whether a set of point in 3-dimensional euclidean space in a cylinder or not
% ==========================================================================================================
% Notice that the set of point need to be a m x 3 matrix, where m is the number of points
% The vectors Cylinder_bot_center and Cylinder_top_center should be two 1 x 3 vectors
% The parameter for the location of cylinder is setted in the function
%                        ' FAME_Parameter_Simple_Cubic '
% 
% ========================================================================================================== 

if Cylinder_radius <= 0
    Point_idx = [];
else
    % Counting the number of points
    % Shifting the points
    Point_bot_set = Point_set - Cylinder_bot_center;
    Point_top_set = Point_set - Cylinder_top_center;
    % Computing the direction of cylinder
    Direct_vector = Cylinder_top_center - Cylinder_bot_center;
    % Normalizing the direction vector
    Direct_vector =       Direct_vector / norm( Direct_vector );
    % Constructing the array for vectorization computing 
%     Direct_vector_array = kron( ones( Point_num , 1 ) , Direct_vector );
    % Computing the distance for each point to the direction vector
%     Point_line_dist   = cross( Point_bot_set , Direct_vector_array ).^2 * ones(3,1);
    Point_line_dist   =  ( Point_bot_set(:,2)*Direct_vector(3) - Point_bot_set(:,3)*Direct_vector(2) ).^2 + ...
                         ( Point_bot_set(:,3)*Direct_vector(1) - Point_bot_set(:,1)*Direct_vector(3) ).^2 + ... 
                         ( Point_bot_set(:,1)*Direct_vector(2) - Point_bot_set(:,2)*Direct_vector(1) ).^2 ;
    % Computing norm of the projection vectors 
    Point_line_bot_proj =    Point_bot_set(:,1)*Direct_vector(1) + Point_bot_set(:,2)*Direct_vector(2) + Point_bot_set(:,3)*Direct_vector(3);
    Point_line_top_proj = - (Point_top_set(:,1)*Direct_vector(1) + Point_top_set(:,2)*Direct_vector(2) + Point_top_set(:,3)*Direct_vector(3) );
%     Point_line_proj     = min( [ Point_line_bot_proj , Point_line_top_proj ],[],2 );
    Point_line_proj     = min( Point_line_bot_proj,Point_line_top_proj );
    % Determining whether Point_set in a cylinder or not
    Point_idx_proj = find( Point_line_proj >= 0 );
    Point_idx_dist = find( Point_line_dist <= Cylinder_radius * Cylinder_radius );
    
    Point_idx = intersect(Point_idx_dist,Point_idx_proj);
end