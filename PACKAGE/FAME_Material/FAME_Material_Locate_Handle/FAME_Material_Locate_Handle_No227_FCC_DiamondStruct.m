function idx = FAME_Material_Locate_Handle_No227_FCC_DiamondStruct( X, Y, Z, a1, a2, a3, sphere_radius, cylinder_radius )
    sphere_center = [0 0 0];
    idx_sphere    = FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center );
    
    cylinder_bot_center{1} = a1;           cylinder_top_center{1} = (a1 + a2 + a3)/3;
    cylinder_bot_center{2} = a2;           cylinder_top_center{2} = (a1 + a2 + a3)/3;
    cylinder_bot_center{3} = a3;           cylinder_top_center{3} = (a1 + a2 + a3)/3;
    cylinder_bot_center{4} = [0;0;0];      cylinder_top_center{4} = (a1 + a2 + a3)/3;
    idx_cylinder = [];
    for i = 1:4
        idx_temp     = FAME_Material_Locate_Cylinder( [X,Y,Z], cylinder_radius, cylinder_bot_center{i}', cylinder_top_center{i}' );
        idx_cylinder = union( idx_cylinder, idx_temp);
    end
    
    idx = union(idx_sphere, idx_cylinder);
end