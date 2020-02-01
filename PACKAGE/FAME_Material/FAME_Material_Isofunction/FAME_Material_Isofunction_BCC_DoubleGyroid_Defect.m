function G = FAME_Material_Isofunction_BCC_DoubleGyroid_Defect( x,y,z,a1,a2,a3,lattice_constant_a, num_airsphere )
Fun_Gyroid_1    = @(x,y,z) sin(2*pi*x/lattice_constant_a).*cos(2*pi*y/lattice_constant_a) + ...
    sin(2*pi*y/lattice_constant_a).*cos(2*pi*z/lattice_constant_a) + ...
    sin(2*pi*z/lattice_constant_a).*cos(2*pi*x/lattice_constant_a);
Fun_Gyroid_2    = @(x,y,z) Fun_Gyroid_1(-x,-y,-z);

G{2} = Fun_Gyroid_1(x,y,z);
G{1} = Fun_Gyroid_2(x,y,z);

%% One air Sphere
% Center1   = 0.375 * a1 + 0.75 * a2 + 0.125 * a3;
Center1   = [1/4 ,-1/8,1/2]';
Fun_Sphere1      = @(x,y,z) sqrt((x - Center1(1)).^2 + (y - Center1(2)).^2 + (z - Center1(3)).^2);
temp1 = -Fun_Sphere1(x,y,z);
G{3} = temp1;
%% Double air Sphere
if num_airsphere == 2
    Center2   = (1 - 0.375) * a1 + (1 - 0.75) * a2 + (1 - 0.125) * a3;
    Fun_Sphere2      = @(x,y,z) sqrt((x - Center2(1)).^2 + (y - Center2(2)).^2 + (z - Center2(3)).^2);
    temp2 = -Fun_Sphere2(x,y,z);
    G{3} = max(temp1, temp2);
end
end