function [sphere_centers,sphere_radius] = ...
            FAME_Plot_Periodic_Primitive_Cell_Sphere(a1,a2,a3,sphere_centers,sphere_radius)
    
    a1 = reshape(a1,3,1);
    a2 = reshape(a2,3,1);
    a3 = reshape(a3,3,1);

    % Determine periodic for sphere

    [ periodic_sphere_centers, periodic_sphere_radius ] = ...
        periodic_determine_sphere(a1,a2,a3,sphere_centers, sphere_radius);
    sphere_centers = [ sphere_centers; periodic_sphere_centers];
%     sphere_radius(end+1:length(periodic_sphere_radius)-1)  = [ sphere_radius   periodic_sphere_radius ];
    sphere_radius(end+1:end+length(periodic_sphere_radius))  =   periodic_sphere_radius ;
end

function [ periodic_sphere_centers, periodic_sphere_radius ] = periodic_determine_sphere(a1,a2,a3,sphere_centers, sphere_radius)

    [ n_sphere,flag ] = input_check(sphere_centers,sphere_radius);

    val = [-1; 0; 1];
    e   = [ 1; 1; 1];
    shift = [kron(e,kron(e,val)) kron(e,kron(val,e)) kron(val,kron(e,e))];

    periodic_sphere_centers = [];
    periodic_sphere_radius  = [];
    for i = 1:n_sphere
        
        for j = 1:size(shift,1)
            sphere_center = sphere_centers(i,:) + shift(j,:);
            point = [a1 a2 a3]*sphere_center';
            d = dist_point2parallelopiped(a1,a2,a3,point);
            if d(1)*d(2) < 1e-3 && d(3)*d(4) < 1e-3 && d(5)*d(6) < 1e-3
                periodic_sphere_centers = [periodic_sphere_centers;sphere_center];
                periodic_sphere_radius  = [periodic_sphere_radius ,sphere_radius(i)];
            end
        end
    end
       
end

function [ n_sphere,flag ] = input_check(sphere_centers,sphere_radius)
    if size(sphere_centers,1) ~= length(sphere_radius)
        error('The number of descriptions for sphere were distinct! Please check the size input data.')
        flag = 1;
    else
        n_sphere = size(sphere_centers,1);
        flag = 0;
    end
end

function d = dist_point2parallelopiped(a1,a2,a3,point)
%     figure(10); hold on
%     x = linspace(-2,2,100);
%     y = linspace(-2,2,100);
%     z = linspace(-2,2,100);
%     [X,Y,Z] = meshgrid(x,y,z);

%     point = [a1 a2 a3]*point';

% for A12
    vec_n = cross(a1,a2);
    fun_test_A12 = @(x,y,z) vec_n(1)*x + vec_n(2)*y + vec_n(3)*z;
%     isosurface(X,Y,Z,fun_test_A12(X,Y,Z),0);

    d(1) = fun_test_A12(point(1),point(2),point(3))/norm(vec_n);
% for A12'
    fun_test_A12p = @(x,y,z) fun_test_A12(x-a3(1),y-a3(2),z-a3(3));
%     isosurface(X,Y,Z,fun_test_A12p(X,Y,Z),0);
    
    d(2) = fun_test_A12p(point(1),point(2),point(3))/norm(vec_n);
% for A13
    vec_n = cross(a1,a3);
    fun_test_A13 = @(x,y,z) vec_n(1)*x + vec_n(2)*y + vec_n(3)*z;
%     isosurface(X,Y,Z,fun_test_A13(X,Y,Z),0);
    
    d(3) = fun_test_A13(point(1),point(2),point(3))/norm(vec_n);
% for A13'
    fun_test_A13p = @(x,y,z) fun_test_A13(x-a2(1),y-a2(2),z-a2(3));
%     isosurface(X,Y,Z,fun_test_A13p(X,Y,Z),0); 
    
    d(4) = fun_test_A13p(point(1),point(2),point(3))/norm(vec_n);
% for A23
    vec_n = cross(a2,a3);
    fun_test_A23 = @(x,y,z) vec_n(1)*x + vec_n(2)*y + vec_n(3)*z;
%     isosurface(X,Y,Z,fun_test_A23(X,Y,Z),0);
    
    d(5) = fun_test_A23(point(1),point(2),point(3))/norm(vec_n);
% for A23'
    fun_test_A23p = @(x,y,z) fun_test_A23(x-a1(1),y-a1(2),z-a1(3));
%     isosurface(X,Y,Z,fun_test_A23p(X,Y,Z),0);   
    
    d(6) = fun_test_A23p(point(1),point(2),point(3))/norm(vec_n);
    
%     camlight left; 
    
    
    
end