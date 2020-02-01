function G = FAME_Material_Isofunction_TR_SGDG_Defect( x,y,z,a1,a2,a3,lattice_constant_a, ncell )
Fun_Gyroid_1    = @(x,y,z) sin(2*pi*x/lattice_constant_a).*cos(2*pi*y/lattice_constant_a) + ...
    sin(2*pi*y/lattice_constant_a).*cos(2*pi*z/lattice_constant_a) + ...
    sin(2*pi*z/lattice_constant_a).*cos(2*pi*x/lattice_constant_a);
Fun_Gyroid_2    = @(x,y,z) Fun_Gyroid_1(-x,-y,-z);

Point_set     = [x,y,z];
coef          = Point_set/([a1,a2,a3]');

G{2} = Fun_Gyroid_1(x,y,z);
G{1} = Fun_Gyroid_2(x,y,z);

n_DG = ncell(2);
n_total = sum(ncell);

G{2}( coef(:,2) > n_DG / n_total ,:) = nan;
G{3} = [];
for i = 1 : n_DG
    Center   = 0.375 * a1 + (0.75 + i - 1) * a2 / n_total  + 0.125 * a3;
% Center   = 0.25 * a1 + (-0.375 + i - 1) * a2 / n_total  + 0.5 * a3;
    Fun_Sphere      = @(x,y,z) sqrt((x - Center(1)).^2 + (y - Center(2)).^2 + (z - Center(3)).^2);
    temp = Fun_Sphere(x,y,z);
    if i == 1
        G{3} = temp;
    else
        G{3} = min(G{3}, temp);
    end
end
G{3} = -G{3};
end