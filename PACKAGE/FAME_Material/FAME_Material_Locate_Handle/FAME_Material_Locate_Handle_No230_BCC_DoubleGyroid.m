function Point_idx = FAME_Material_Locate_Handle_No230_BCC_DoubleGyroid( X, Y, Z, a1, a2, a3, iso_num, lattice_constant )
Fun_Gyroid1    = @(x,y,z) sin(2*pi*x/lattice_constant.a).*cos(2*pi*y/lattice_constant.a) + ...
                          sin(2*pi*y/lattice_constant.a).*cos(2*pi*z/lattice_constant.a) + ...
                          sin(2*pi*z/lattice_constant.a).*cos(2*pi*x/lattice_constant.a);
Fun_Gyroid2    = @(x,y,z) Fun_Gyroid1(-x,-y,-z);

Point_idx = union( find( Fun_Gyroid1(X,Y,Z)  >  iso_num), find( Fun_Gyroid2(X,Y,Z)  >  iso_num) );      
end