function G = FAME_Material_Isofunction_TR_SGDG( x,y,z,a1,a2,a3,lattice_constant_a, ncell)
    Fun_Gyroid_1    = @(x,y,z) sin(2*pi*x/lattice_constant_a).*cos(2*pi*y/lattice_constant_a) + ...
                               sin(2*pi*y/lattice_constant_a).*cos(2*pi*z/lattice_constant_a) + ...
                               sin(2*pi*z/lattice_constant_a).*cos(2*pi*x/lattice_constant_a);
    Fun_Gyroid_2    = @(x,y,z) Fun_Gyroid_1(-x,-y,-z);
    
    Point_set     = [x,y,z];
    coef          = Point_set/([a1,a2,a3]');
    
    G{1} = Fun_Gyroid_1(x,y,z);
    G{2} = Fun_Gyroid_2(x,y,z);
    
    n_DG = ncell(2);
    n_total = sum(ncell);
    G{2}( coef(:,2) > n_DG / n_total ,:) = nan;
end



