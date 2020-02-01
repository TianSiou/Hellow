function idx = FAME_Material_Locate_Handle_No215_SC_Fe4C( X, Y, Z, a1, a2, a3, sphere_radius, cylinder_radius )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% When x2 = 1/4, the iron atoms are at the positions of the face-centered cubic lattice.
% In Fe4C, x2 is about 0.265. (that is, the lattice constant a is about 4*0.265)
% See detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/AB4_cP5_215_a_e.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x2 = 0.265;
    % Given sphere centers
    sphere_center_C     = [0; 0; 0];
    sphere_center_Fe_1  =     x2*a1 +     x2*a2 +     x2*a3;
    sphere_center_Fe_2  = (1-x2)*a1 + (1-x2)*a2 +     x2*a3;
    sphere_center_Fe_3  = (1-x2)*a1 +     x2*a2 + (1-x2)*a3;
    sphere_center_Fe_4  =     x2*a1 + (1-x2)*a2 + (1-x2)*a3;
    % Determine indices inner material
    idx_sphere_Fe_1  =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_Fe_1' );
    idx_sphere_Fe_2  =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_Fe_2' );
    idx_sphere_Fe_3  =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_Fe_3' );
    idx_sphere_Fe_4  =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_Fe_4' );
    
    idx_sphere_C   =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_C' );
    idx_sphere_Fe  =  union(idx_sphere_Fe_1,union(idx_sphere_Fe_2,union(idx_sphere_Fe_3,idx_sphere_Fe_4))) ;
    
    % Given top centers and bot centers of cylinders
    cylinder_bot_center_C_1 = [0; 0; 0];   cylinder_top_center_C_1 =            sphere_center_Fe_1 /2;
    cylinder_bot_center_C_2 = a1 + a2;     cylinder_top_center_C_2 = (a1 + a2 + sphere_center_Fe_2)/2;
    cylinder_bot_center_C_3 = a1 + a3;     cylinder_top_center_C_3 = (a1 + a3 + sphere_center_Fe_3)/2;
    cylinder_bot_center_C_4 = a2 + a3;     cylinder_top_center_C_4 = (a2 + a3 + sphere_center_Fe_4)/2;
    cylinder_bot_center_Fe_1 = sphere_center_Fe_1;   cylinder_top_center_Fe_1 =            sphere_center_Fe_1 /2;
    cylinder_bot_center_Fe_2 = sphere_center_Fe_2;   cylinder_top_center_Fe_2 = (a1 + a2 + sphere_center_Fe_2)/2;
    cylinder_bot_center_Fe_3 = sphere_center_Fe_3;   cylinder_top_center_Fe_3 = (a1 + a3 + sphere_center_Fe_3)/2;
    cylinder_bot_center_Fe_4 = sphere_center_Fe_4;   cylinder_top_center_Fe_4 = (a2 + a3 + sphere_center_Fe_4)/2;
    % Determine indices inner material
    idx_cylinder_C_1   =  FAME_Material_Locate_Cylinder( [X,Y,Z], cylinder_radius, cylinder_bot_center_C_1', cylinder_top_center_C_1' );
    idx_cylinder_C_2   =  FAME_Material_Locate_Cylinder( [X,Y,Z], cylinder_radius, cylinder_bot_center_C_2', cylinder_top_center_C_2' );
    idx_cylinder_C_3   =  FAME_Material_Locate_Cylinder( [X,Y,Z], cylinder_radius, cylinder_bot_center_C_3', cylinder_top_center_C_3' );
    idx_cylinder_C_4   =  FAME_Material_Locate_Cylinder( [X,Y,Z], cylinder_radius, cylinder_bot_center_C_4', cylinder_top_center_C_4' );
    idx_cylinder_Fe_1  =  FAME_Material_Locate_Cylinder( [X,Y,Z], cylinder_radius, cylinder_bot_center_Fe_1', cylinder_top_center_Fe_1' );
    idx_cylinder_Fe_2  =  FAME_Material_Locate_Cylinder( [X,Y,Z], cylinder_radius, cylinder_bot_center_Fe_2', cylinder_top_center_Fe_2' );
    idx_cylinder_Fe_3  =  FAME_Material_Locate_Cylinder( [X,Y,Z], cylinder_radius, cylinder_bot_center_Fe_3', cylinder_top_center_Fe_3' );
    idx_cylinder_Fe_4  =  FAME_Material_Locate_Cylinder( [X,Y,Z], cylinder_radius, cylinder_bot_center_Fe_4', cylinder_top_center_Fe_4' );
    
    idx_cylinder_C   =  union(  idx_cylinder_C_1,  union(idx_cylinder_C_2,  union(idx_cylinder_C_3,  idx_cylinder_C_4)) );
    idx_cylinder_Fe  =  union( idx_cylinder_Fe_1, union(idx_cylinder_Fe_2, union(idx_cylinder_Fe_3, idx_cylinder_Fe_4)) );
    
    % Merge indices
    idx{1}  =  union( idx_sphere_C,  idx_cylinder_C); % indices inner C
    idx{2}  =  union(idx_sphere_Fe, idx_cylinder_Fe); % indices inner Fe
end