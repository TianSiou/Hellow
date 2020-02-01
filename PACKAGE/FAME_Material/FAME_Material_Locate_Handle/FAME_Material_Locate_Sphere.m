function [ Point_idx ] = FAME_Material_Locate_Sphere( Point_set, Sphere_radius, Sphere_center )
%% This function is used to determine whether a set of point in 3-dimensional euclidean space in seven spheres or not
% ==========================================================================================================
% Notice that the set of point need to be a m x 3 matrix, where m is the number of points
% The vector Sphere_center should be a 1 x 3 vector.  
% The parameter for the location of sphere is setting in the function
%                        ' FAME_Material_Parameter_Simple_Cube '
% 
% ==========================================================================================================
% Sum = ( Point_set - Sphere_center ).^2 *ones(3,1);
Sum = ( Point_set(:,1) - Sphere_center(1) ).^2 + ( Point_set(:,2) - Sphere_center(2) ).^2 + ( Point_set(:,3) - Sphere_center(3) ).^2;
Point_idx =  find(Sum <= Sphere_radius*Sphere_radius);
end