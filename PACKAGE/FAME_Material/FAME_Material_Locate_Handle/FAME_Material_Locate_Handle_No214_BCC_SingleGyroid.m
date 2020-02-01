function Point_idx = FAME_Material_Locate_Handle_No214_BCC_SingleGyroid( X, Y, Z, a1, a2, a3, iso_num, lattice_constant )
Fun_Gyroid1    = @(x,y,z) sin(2*pi*x/lattice_constant).*cos(2*pi*y/lattice_constant) + ...
                          sin(2*pi*y/lattice_constant).*cos(2*pi*z/lattice_constant) + ...
                          sin(2*pi*z/lattice_constant).*cos(2*pi*x/lattice_constant);

Point_idx =  find( Fun_Gyroid1(X,Y,Z)  >  iso_num);      
end