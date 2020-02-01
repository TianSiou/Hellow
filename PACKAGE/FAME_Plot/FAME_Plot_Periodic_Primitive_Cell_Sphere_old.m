function [sphere_centers,sphere_radius] = ...
            FAME_Plot_Periodic_Primitive_Cell_Sphere(a1,a2,a3,sphere_centers,sphere_radius)
    
    a1 = reshape(a1,3,1);
    a2 = reshape(a2,3,1);
    a3 = reshape(a3,3,1);

    % Determine periodic for sphere

    [ periodic_sphere_centers, periodic_sphere_radius ] = ...
        periodic_determine_sphere(a1,a2,a3,sphere_centers, sphere_radius);
    sphere_centers = [ sphere_centers; periodic_sphere_centers];
    sphere_radius  = [ sphere_radius   periodic_sphere_radius ];

end

function [ periodic_sphere_centers, periodic_sphere_radius ] = periodic_determine_sphere(a1,a2,a3,sphere_centers, sphere_radius)

    [ n_sphere,flag ] = input_check(sphere_centers,sphere_radius);

    aa1 = [1; 0; 0];
    aa2 = [0; 1; 0];
    aa3 = [0; 0; 1];

    periodic_sphere_centers = [];
    periodic_sphere_radius  = [];
    for i = 1:n_sphere
        sphere_center = sphere_centers(i,:);
        
        dist = dist_point2parallelopiped(a1,a2,a3,sphere_center);
        idx  = find(dist < sphere_radius(i) + 1e-6);

    %% one toutch
        if length(idx) == 1
            shift = [aa3';-aa3';aa2';-aa2';aa1';-aa1'];
            for j = 1:6             
                if ismember(j,idx)
                    periodic_sphere_centers = [ periodic_sphere_centers; sphere_center+shift(j,:) ];
                    periodic_sphere_radius = [ periodic_sphere_radius sphere_radius(i) ];
                end
            end
        end
    %% two toutch
        if length(idx) == 2
            shift = { aa3'  aa2'  aa3'+aa2';
                      aa3' -aa2'  aa3'-aa2';
                      aa3'  aa1'  aa3'+aa1';
                      aa3' -aa1'  aa3'-aa1';
                     -aa3'  aa2' -aa3'+aa2';
                     -aa3' -aa2' -aa3'-aa2';
                     -aa3'  aa1' -aa3'+aa1';
                     -aa3' -aa1' -aa3'-aa1';
                      aa2'  aa1'  aa2'+aa1';
                      aa2' -aa1'  aa2'-aa1';
                     -aa2'  aa1' -aa2'+aa1';
                     -aa2' -aa1' -aa2'-aa1'};     
            pair = {[1 3],[1 4],[1 5],[1 6],[2 3],[2 4],[2 5],[2 6],[3 5],[3 6],[4 5],[4 6]};
            for j = 1:length(pair)
                if ismember(pair{j},idx)
                    periodic_sphere_centers = [ periodic_sphere_centers; sphere_center+shift{j,1}; sphere_center+shift{j,2}; sphere_center+shift{j,3} ];
                    periodic_sphere_radius = [ periodic_sphere_radius sphere_radius(i) sphere_radius(i) sphere_radius(i) ];
                end
            end
        end
    %% three toutch
        if length(idx) == 3
            shift = { aa3'  aa2'  aa1'  aa3'+aa2'  aa3'+aa1'  aa2'+aa1'  aa3'+aa2'+aa1'; % 1 3 5
                      aa3'  aa2' -aa1'  aa3'+aa2'  aa3'-aa1'  aa2'-aa1'  aa3'+aa2'-aa1'; % 1 3 6
                      aa3' -aa2'  aa1'  aa3'-aa2'  aa3'+aa1' -aa2'+aa1'  aa3'-aa2'+aa1'; % 1 4 5
                     -aa3'  aa2'  aa1' -aa3'+aa2' -aa3'+aa1'  aa2'+aa1' -aa3'+aa2'+aa1'; % 2 3 5
                     -aa3' -aa2'  aa1' -aa3'-aa2' -aa3'+aa1' -aa2'+aa1' -aa3'-aa2'+aa1'; % 2 4 5
                     -aa3'  aa2' -aa1' -aa3'+aa2'  aa3'-aa1'  aa2'-aa1' -aa3'+aa2'-aa1'; % 2 3 6
                      aa3' -aa2' -aa1'  aa3'-aa2'  aa3'-aa1' -aa2'-aa1'  aa3'-aa2'-aa1'; % 1 4 6
                     -aa3' -aa2' -aa1' -aa3'-aa2' -aa3'-aa1' -aa2'-aa1' -aa3'-aa2'-aa1'};% 2 4 6
            pair = {[1 3 5],[1 3 6],[1 4 5],[2 3 5],[2 4 5],[2 3 6],[1 4 6],[2 4 6]};
            for j = 1:length(pair)
                if ismember(pair{j},idx)
                    periodic_sphere_centers = [ periodic_sphere_centers; sphere_center+shift{j,1}; sphere_center+shift{j,2}; sphere_center+shift{j,3}; sphere_center+shift{j,4};
                                                                         sphere_center+shift{j,5}; sphere_center+shift{j,6}; sphere_center+shift{j,7} ];
                    periodic_sphere_radius = [ periodic_sphere_radius sphere_radius(i) sphere_radius(i) sphere_radius(i) sphere_radius(i) sphere_radius(i) sphere_radius(i) sphere_radius(i) ];
                end
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

    point = [a1 a2 a3]*point';

% for A12
    vec_n = cross(a1,a2);
    fun_test_A12 = @(x,y,z) vec_n(1)*x + vec_n(2)*y + vec_n(3)*z;
%     isosurface(X,Y,Z,fun_test_A12(X,Y,Z),0);

    d(1) = abs(fun_test_A12(point(1),point(2),point(3)))/norm(vec_n);
% for A12'
    fun_test_A12p = @(x,y,z) fun_test_A12(x-a3(1),y-a3(2),z-a3(3));
%     isosurface(X,Y,Z,fun_test_A12p(X,Y,Z),0);
    
    d(2) = abs(fun_test_A12p(point(1),point(2),point(3)))/norm(vec_n);
% for A13
    vec_n = cross(a1,a3);
    fun_test_A13 = @(x,y,z) vec_n(1)*x + vec_n(2)*y + vec_n(3)*z;
%     isosurface(X,Y,Z,fun_test_A13(X,Y,Z),0);
    
    d(3) = abs(fun_test_A13(point(1),point(2),point(3)))/norm(vec_n);
% for A13'
    fun_test_A13p = @(x,y,z) fun_test_A13(x-a2(1),y-a2(2),z-a2(3));
%     isosurface(X,Y,Z,fun_test_A13p(X,Y,Z),0); 
    
    d(4) = abs(fun_test_A13p(point(1),point(2),point(3)))/norm(vec_n);
% for A23
    vec_n = cross(a2,a3);
    fun_test_A23 = @(x,y,z) vec_n(1)*x + vec_n(2)*y + vec_n(3)*z;
%     isosurface(X,Y,Z,fun_test_A23(X,Y,Z),0);
    
    d(5) = abs(fun_test_A23(point(1),point(2),point(3)))/norm(vec_n);
% for A23'
    fun_test_A23p = @(x,y,z) fun_test_A23(x-a1(1),y-a1(2),z-a1(3));
%     isosurface(X,Y,Z,fun_test_A23p(X,Y,Z),0);   
    
    d(6) = abs(fun_test_A23p(point(1),point(2),point(3)))/norm(vec_n);
    
%     camlight left; 
    
    
    
end