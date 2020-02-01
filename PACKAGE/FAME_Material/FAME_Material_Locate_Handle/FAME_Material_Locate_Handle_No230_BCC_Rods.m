function idx = FAME_Material_Locate_Handle_No230_BCC_Rods( X, Y, Z, a1, a2, a3, cylinder_radius )    
    b1 = a2 + a3;
    b2 = a3 + a1;
    b3 = a1 + a2;

    cylinder_bot_center_r  = [0; 0; 0];         cylinder_top_center_r  = b1 + b2 + b3;
    cylinder_bot_center_g1 = b1 + b3/2;         cylinder_top_center_g1 = b1/2 + b2/2 + b3;
    cylinder_bot_center_g2 = b1/2 + b2/2;       cylinder_top_center_g2 = b2 + b3/2;
    cylinder_bot_center_y1 = b1 + b2/2;         cylinder_top_center_y1 = b1/2 + b3/2;
    cylinder_bot_center_y2 = b1/2 + b2 + b3/2;  cylinder_top_center_y2 = b2/2 + b3;
    cylinder_bot_center_b1 = b1 + b2/2 + b3/2;  cylinder_top_center_b1 = b1/2 + b2;
    cylinder_bot_center_b2 = b1/2 + b3;         cylinder_top_center_b2 = b2/2 + b3/2;
    
    idx_cylinder_r = FAME_Material_Locate_Cylinder( [X,Y,Z], cylinder_radius, cylinder_bot_center_r', cylinder_top_center_r' );
    idx_cylinder_g = union( FAME_Material_Locate_Cylinder( [X,Y,Z], cylinder_radius, cylinder_bot_center_g1', cylinder_top_center_g1' ),...
                            FAME_Material_Locate_Cylinder( [X,Y,Z], cylinder_radius, cylinder_bot_center_g2', cylinder_top_center_g2' ) );
    idx_cylinder_y = union( FAME_Material_Locate_Cylinder( [X,Y,Z], cylinder_radius, cylinder_bot_center_y1', cylinder_top_center_y1' ),...
                            FAME_Material_Locate_Cylinder( [X,Y,Z], cylinder_radius, cylinder_bot_center_y2', cylinder_top_center_y2' ) );
    idx_cylinder_b = union( FAME_Material_Locate_Cylinder( [X,Y,Z], cylinder_radius, cylinder_bot_center_b1', cylinder_top_center_b1' ),...
                            FAME_Material_Locate_Cylinder( [X,Y,Z], cylinder_radius, cylinder_bot_center_b2', cylinder_top_center_b2' ) );
    
    idx{1} = idx_cylinder_r;
    idx{2} = idx_cylinder_g;
    idx{3} = idx_cylinder_y;
    idx{4} = idx_cylinder_b;
end