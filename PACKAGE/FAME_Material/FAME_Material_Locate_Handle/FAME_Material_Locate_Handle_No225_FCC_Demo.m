function idx = FAME_Material_Locate_Handle_No225_FCC_Demo( X, Y, Z, a1, a2, a3, sphere_radius, cylinder_radius )
    sphere_center = (a1 + a2 + a3)/2;
    idx_sphere    = FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center );
    
    cylinder_bot_center_1 = (a2 + a3)/2; cylinder_top_center_1 = (a2 + a3)/2 + a1;
    cylinder_bot_center_2 = (a3 + a1)/2; cylinder_top_center_2 = (a3 + a1)/2 + a2;
    cylinder_bot_center_3 = (a1 + a2)/2; cylinder_top_center_3 = (a1 + a2)/2 + a3;
    idx_cylinder_1        = FAME_Material_Locate_Cylinder( [X,Y,Z], cylinder_radius, cylinder_bot_center_1', cylinder_top_center_1' );
    idx_cylinder_2        = FAME_Material_Locate_Cylinder( [X,Y,Z], cylinder_radius, cylinder_bot_center_2', cylinder_top_center_2' );
    idx_cylinder_3        = FAME_Material_Locate_Cylinder( [X,Y,Z], cylinder_radius, cylinder_bot_center_3', cylinder_top_center_3' );
    idx_cylinder          = union( idx_cylinder_1, union(idx_cylinder_2,idx_cylinder_3) );
    
    idx = union(idx_sphere, idx_cylinder);
end