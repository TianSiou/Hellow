function G = FAME_Material_Isofunction_BCC_DoubleGyroid( x,y,z,a1,a2,a3,lattice_constant_a )
    Fun_Gyroid_1    = @(x,y,z) sin(2*pi*x/lattice_constant_a).*cos(2*pi*y/lattice_constant_a) + ...
                               sin(2*pi*y/lattice_constant_a).*cos(2*pi*z/lattice_constant_a) + ...
                               sin(2*pi*z/lattice_constant_a).*cos(2*pi*x/lattice_constant_a);
    Fun_Gyroid_2    = @(x,y,z) Fun_Gyroid_1(-x,-y,-z);
    
    G{1} = Fun_Gyroid_1(x,y,z);
    G{2} = Fun_Gyroid_2(x,y,z);
end