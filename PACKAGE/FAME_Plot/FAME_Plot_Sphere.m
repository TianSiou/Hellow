function h = FAME_Plot_Sphere(a1,a2,a3,center, radius)
%     x_min = min([a1(1), a2(1), a3(1)]);  x_max = max([a1(1), a2(1), a3(1)]);
%     y_min = min([a1(2), a2(3), a3(2)]);  y_max = max([a1(2), a2(2), a3(2)]);
%     z_min = min([a1(3), a2(3), a3(3)]);  z_max = max([a1(3), a2(3), a3(3)]);
%     
%     n = 60;
%     x = linspace(x_min, x_max, n);
%     y = linspace(y_min, y_max, n);
%     z = linspace(z_min, z_max, n);
%     [X,Y,Z] = meshgrid(x,y,z);

    center = center*[a1';a2';a3'];
%     fun = @(x,y,z,c) ( x-c(1) )^2 + ( y-c(2) )^2 + ( z-c(3) )^2 ;
    for i = 1:length(radius)
        [X_sp,Y_sp,Z_sp] = sphere(20);
        h(i) = surf(radius(i)*X_sp+center(i,1), radius(i)*Y_sp+center(i,2), radius(i)*Z_sp+center(i,3) );
%         val = fun(X,Y,Z,center(i,:));
%         h(i) = isosurface(X,Y,Z,val,radius(i)^2);
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

    d(1) = abs(fun_test_A12(point(1,:),point(2,:),point(3,:)))/norm(vec_n);
% for A12'
    fun_test_A12p = @(x,y,z) fun_test_A12(x-a3(1),y-a3(2),z-a3(3));
%     isosurface(X,Y,Z,fun_test_A12p(X,Y,Z),0);
    
    d(2) = fun_test_A12p(point(1,:),point(2,:),point(3,:));
% for A13
    vec_n = cross(a1,a3);
    fun_test_A13 = @(x,y,z) vec_n(1)*x + vec_n(2)*y + vec_n(3)*z;
%     isosurface(X,Y,Z,fun_test_A13(X,Y,Z),0);
    
    d(3) = fun_test_A13(point(1,:),point(2,:),point(3,:));
% for A13'
    fun_test_A13p = @(x,y,z) fun_test_A13(x-a2(1),y-a2(2),z-a2(3));
%     isosurface(X,Y,Z,fun_test_A13p(X,Y,Z),0); 
    
    d(4) = fun_test_A13p(point(1,:),point(2,:),point(3,:));
% for A23
    vec_n = cross(a2,a3);
    fun_test_A23 = @(x,y,z) vec_n(1)*x + vec_n(2)*y + vec_n(3)*z;
    
    d(5) = fun_test_A23(point(1,:),point(2,:),point(3,:));
% for A23'
    fun_test_A23p = @(x,y,z) fun_test_A23(x-a1(1),y-a1(2),z-a1(3));
    
    d(6) = fun_test_A23p(point(1,:),point(2,:),point(3,:)); 
end