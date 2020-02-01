function [cylinder_bot_centers,cylinder_top_centers,cylinder_radius] = ...
            FAME_Plot_Periodic_Primitive_Cell_Cylinder(a1,a2,a3,cylinder_bot_centers, cylinder_top_centers, cylinder_radius)
% size(cylinder_bot_centers,1)
% size(cylinder_top_centers,1)
% length(cylinder_radius)
    a1 = reshape(a1,3,1);
    a2 = reshape(a2,3,1);
    a3 = reshape(a3,3,1);
    
    % Determine periodic for cylinder
    [ periodic_cylinder_bot_centers, periodic_cylinder_top_centers, periodic_cylinder_radius ] = ...
        periodic_determine_cylinder(a1,a2,a3,cylinder_bot_centers, cylinder_top_centers, cylinder_radius);
    cylinder_bot_centers = [ cylinder_bot_centers; periodic_cylinder_bot_centers];
    cylinder_top_centers = [ cylinder_top_centers; periodic_cylinder_top_centers];
    cylinder_radius      = [ cylinder_radius       periodic_cylinder_radius     ];

end

function [ periodic_cylinder_bot_centers, periodic_cylinder_top_centers, periodic_cylinder_radius ] = periodic_determine_cylinder(a1,a2,a3,cylinder_bot_centers, cylinder_top_centers,cylinder_radius)

    [ n_cylinder, flag ] = input_check(cylinder_bot_centers, cylinder_top_centers,cylinder_radius);
    
    tol = 1e-2;
    
    aa1 = [1; 0; 0];
    aa2 = [0; 1; 0];
    aa3 = [0; 0; 1];

    periodic_cylinder_bot_centers = [];
    periodic_cylinder_top_centers = [];
    periodic_cylinder_radius      = [];
    for i = 1:n_cylinder
        cylinder_bot_center = cylinder_bot_centers(i,:);
        cylinder_top_center = cylinder_top_centers(i,:);

        direction = [a1 a2 a3]*(cylinder_top_center - cylinder_bot_center)';
        parallel_vector  = [        norm(cross(direction,a1)),        norm(cross(direction,a2)),        norm(cross(direction,a3)) ];
        parallel_plane   = [ abs(dot(direction,cross(a1,a2))), abs(dot(direction,cross(a1,a3))), abs(dot(direction,cross(a2,a3))) ];

        case_num_I  = find( parallel_vector < 1e-6 ); % 1:  vec//a1; 2:  vec//a2; 3:  vec//a3; 
        case_num_II = find( parallel_plane  < 1e-6 ); % 1: vec//A12; 2: vec//A13; 3: vec//A23; 

        shift = [];

        if isempty(case_num_I) == 0
            point = (cylinder_bot_center + cylinder_top_center)/2;
            switch case_num_I
                case 1
                    shift_temp = {aa3', -aa3'; aa2', -aa2'};
                    dist_12 = dist_point2plane(a1,a2,a3,point,'12');
                    dist_13 = dist_point2plane(a1,a2,a3,point,'13');
                    temp_12 = find( dist_12 < cylinder_radius(i)+tol );
                    temp_13 = find( dist_13 < cylinder_radius(i)+tol );
                    if length(temp_12) + length(temp_13) == 1
                        shift   = [shift; shift_temp{1,temp_12}; shift_temp{2,temp_13}];
                    elseif length(temp_12) + length(temp_13) == 2
                        shift   = [shift; shift_temp{1,temp_12}; shift_temp{2,temp_13}; shift_temp{1,temp_12}+shift_temp{2,temp_13}];
                    end

                case 2
                    shift_temp = {aa3', -aa3'; aa1', -aa1'};
                    dist_12 = dist_point2plane(a1,a2,a3,point,'12');
                    dist_23 = dist_point2plane(a1,a2,a3,point,'23');
                    temp_12 = find( dist_12 < cylinder_radius(i)+tol );
                    temp_23 = find( dist_23 < cylinder_radius(i)+tol );
                    if length(temp_12) + length(temp_23) == 1
                        shift   = [shift; shift_temp{1,temp_12}; shift_temp{2,temp_23}];
                    elseif length(temp_12) + length(temp_23) == 2
                        shift   = [shift; shift_temp{1,temp_12}; shift_temp{2,temp_23}; shift_temp{1,temp_12}+shift_temp{2,temp_23}];
                    end

                case 3
                    shift_temp = {aa2', -aa2'; aa1', -aa1'};
                    dist_13 = dist_point2plane(a1,a2,a3,point,'13');
                    dist_23 = dist_point2plane(a1,a2,a3,point,'23');
                    temp_13 = find( dist_13 < cylinder_radius(i)+tol );
                    temp_23 = find( dist_23 < cylinder_radius(i)+tol );
                    if length(temp_13) + length(temp_23) == 1
                        shift   = [shift; shift_temp{1,temp_13}; shift_temp{2,temp_23}];
                    elseif length(temp_13) + length(temp_23) == 2
                        shift   = [shift; shift_temp{1,temp_13}; shift_temp{2,temp_23}; shift_temp{1,temp_13}+shift_temp{2,temp_23}];
                    end

            end

        elseif isempty(case_num_II) == 0
            point = [a1 a2 a3]*(cylinder_bot_center + cylinder_top_center)'/2;
            switch case_num_II
                case 1
                    shift_temp = {aa3', -aa3'};
                    dist_12 = dist_point2plane(a1,a2,a3,point,'12');
                    temp_12 = find( dist_12 < cylinder_radius(i)+tol );
                    shift   = [shift; shift_temp{1,temp_12}];

                case 2
                    shift_temp = {aa2', -aa2'};
                    dist_13 = dist_point2plane(a1,a2,a3,point,'13');
                    temp_13 = find( dist_13 < cylinder_radius(i)+tol );
                    shift   = [shift; shift_temp{1,temp_13}];

                case 3
                    shift_temp = {aa1', -aa1'};
                    dist_23 = dist_point2plane(a1,a2,a3,point,'23');
                    temp_23 = find( dist_23 < cylinder_radius(i)+tol );
                    shift   = [shift; shift_temp{1,temp_23}];

            end
        end

        for j = 1:size(shift,1)
            periodic_cylinder_bot_centers = [ periodic_cylinder_bot_centers; cylinder_bot_center+shift(j,:) ];
            periodic_cylinder_top_centers = [ periodic_cylinder_top_centers; cylinder_top_center+shift(j,:) ];
            periodic_cylinder_radius      = [ periodic_cylinder_radius cylinder_radius(i) ];
        end
    end
end

function [ n_cylinder,flag ] = input_check(cylinder_bot_centers, cylinder_top_centers,cylinder_radius)

    if size(cylinder_bot_centers,1) ~= size(cylinder_top_centers,1) || size(cylinder_bot_centers,1) ~= length(cylinder_radius) || size(cylinder_top_centers,1) ~= length(cylinder_radius)
        error('The number of descriptions for cylinder were distinct! Please check the size input data.')
        flag = 1;
    else
        n_cylinder = size(cylinder_bot_centers,1);
        flag = 0;
    end
end

function d = dist_point2plane(a1,a2,a3,point,plane_num)

switch plane_num
    case '12'
% for A12
    vec_n = cross(a1,a2);
    fun_test_A12  = @(x,y,z) vec_n(1)*x + vec_n(2)*y + vec_n(3)*z;
    fun_test_A12p = @(x,y,z) fun_test_A12(x-a3(1),y-a3(2),z-a3(3));
    
    d = [ abs(fun_test_A12(point(1),point(2),point(3)))/norm(vec_n), abs(fun_test_A12p(point(1),point(2),point(3)))/norm(vec_n) ];
    case '13'
% for A13
    vec_n = cross(a1,a3);
    fun_test_A13  = @(x,y,z) vec_n(1)*x + vec_n(2)*y + vec_n(3)*z;
    fun_test_A13p = @(x,y,z) fun_test_A13(x-a2(1),y-a2(2),z-a2(3));
    
    d = [ abs(fun_test_A13(point(1),point(2),point(3)))/norm(vec_n), abs(fun_test_A13p(point(1),point(2),point(3)))/norm(vec_n) ];
    case '23'
% for A23
    vec_n = cross(a2,a3);
    fun_test_A23  = @(x,y,z) vec_n(1)*x + vec_n(2)*y + vec_n(3)*z;
    fun_test_A23p = @(x,y,z) fun_test_A23(x-a1(1),y-a1(2),z-a1(3));
    
    d = [ abs(fun_test_A23(point(1),point(2),point(3)))/norm(vec_n), abs(fun_test_A23p(point(1),point(2),point(3)))/norm(vec_n) ];
end
    
    
end