function Point_idx = FAME_Material_Locate_Handle_User_Defined_Isosurface( X, Y, Z, a1, a2, a3, iso_num, fun_handle )
    
    Point_idx = [];

    for i = 1:length(fun_handle)
        temp_fun_handle = fun_handle{i};
        Point_idx = union( Point_idx, find( temp_fun_handle(X,Y,Z,a1,a2,a3) > iso_num(i)) );
    end
    
end