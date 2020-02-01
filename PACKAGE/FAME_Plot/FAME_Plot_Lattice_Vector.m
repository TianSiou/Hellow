function FAME_Plot_Lattice_Vector(lattice_constant)
    a1 = [ lattice_constant.a; 0; 0 ];
    a2 = [ lattice_constant.b*cos(lattice_constant.theta1); lattice_constant.b*sin(lattice_constant.theta1); 0];
    a3 = [ lattice_constant.c*cos(lattice_constant.theta2); 
           lattice_constant.c*( cos(lattice_constant.theta3) - cos(lattice_constant.theta1)*cos(lattice_constant.theta2) )/sin(lattice_constant.theta1);
           lattice_constant.c*sqrt( 1 - cos(lattice_constant.theta1)^2 - cos(lattice_constant.theta2)^2 - cos(lattice_constant.theta3)^2 + ...
                      2*cos(lattice_constant.theta1)*cos(lattice_constant.theta2)*cos(lattice_constant.theta3) )/sin(lattice_constant.theta1)];
    figure(1); hold on
    
    plot_Parallelepiped(a1,a2,a3)
end