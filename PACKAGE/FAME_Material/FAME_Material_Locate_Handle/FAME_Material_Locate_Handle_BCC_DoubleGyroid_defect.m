function Point_idx = FAME_Material_Locate_Handle_BCC_DoubleGyroid_defect( X, Y, Z, a1, a2, a3, iso_num, lattice_constant )
Fun_Gyroid1    = @(x,y,z) sin(2*pi*x/lattice_constant.a).*cos(2*pi*y/lattice_constant.a) + ...
                          sin(2*pi*y/lattice_constant.a).*cos(2*pi*z/lattice_constant.a) + ...
                          sin(2*pi*z/lattice_constant.a).*cos(2*pi*x/lattice_constant.a);
Fun_Gyroid2    = @(x,y,z) Fun_Gyroid1(-x,-y,-z);

sphere_center = [0.25 , -0.125 , 0.5] ;
Fun_shpere   = @(x,y,z) (x-sphere_center(1)).^2+(y-sphere_center(2)).^2+(z-sphere_center(3)).^2  ; 
Point_idx_sphere = find(Fun_shpere(X,Y,Z)<0.15^2);

Point_idx =  setdiff( union( find( Fun_Gyroid1(X,Y,Z)  >  iso_num), find( Fun_Gyroid2(X,Y,Z)  >  iso_num) ), Point_idx_sphere);      
end