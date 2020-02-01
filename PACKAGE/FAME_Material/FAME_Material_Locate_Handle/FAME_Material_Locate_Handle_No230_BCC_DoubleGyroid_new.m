function Fun_Gyroid = FAME_Material_Locate_Handle_No230_BCC_DoubleGyroid_new( lattice_constant )
    
    Fun_Gyroid{1}    = @(x,y,z,a1,a2,a3) sin(2*pi*x/lattice_constant.a).*cos(2*pi*y/lattice_constant.a) + ...
                                         sin(2*pi*y/lattice_constant.a).*cos(2*pi*z/lattice_constant.a) + ...
                                         sin(2*pi*z/lattice_constant.a).*cos(2*pi*x/lattice_constant.a);
    Fun_Gyroid{2}    = @(x,y,z,a1,a2,a3) Fun_Gyroid{1}(-x,-y,-z);
end