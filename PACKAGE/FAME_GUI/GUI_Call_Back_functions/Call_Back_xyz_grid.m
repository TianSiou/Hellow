function Call_Back_xyz_grid(hobj,event,h)

    x_grid_num = get( h.edit_domain.grid_num(1), 'string' );
    y_grid_num = get( h.edit_domain.grid_num(2), 'string' );
    z_grid_num = get( h.edit_domain.grid_num(3), 'string' );
    
    [Popt] = FAME_GUI_Info_Get(h);
    
    x_mesh_len = num2str( Popt.mesh.x_edge_len / str2num(x_grid_num) );
    y_mesh_len = num2str( Popt.mesh.y_edge_len / str2num(y_grid_num) );
    z_mesh_len = num2str( Popt.mesh.z_edge_len / str2num(z_grid_num) );
    
    FAME_GUI_Info_Set('grid_numbers', {x_grid_num, y_grid_num, z_grid_num}, h)
    FAME_GUI_Info_Set('mesh_length' , {x_mesh_len, y_mesh_len, z_mesh_len}, h)
    
end