function Call_Back_shape_parameter(hobj,event,h)
    
    hobj_idx   = get(hobj,'value');
    hobj_value = str2num(get(hobj,'string'));

    [Popt] = FAME_GUI_Info_Get(h);

    if     hobj_idx == 1
        sphere_radius   = hobj_value;
        cylinder_radius = Popt.material.cylinder_radius;
    elseif hobj_idx == 2
        sphere_radius   = Popt.material.sphere_radius;
        cylinder_radius = hobj_value;
    end
%% show 3D 對新的 edge length 做圖
    lattice_vec_a = [ reshape(Popt.lattice.lattice_vector_a1,3,1), reshape(Popt.lattice.lattice_vector_a2,3,1), reshape(Popt.lattice.lattice_vector_a3,3,1)];

    FAME_Main_User_Defined_gui(Popt.material.file_name,sphere_radius,cylinder_radius, ...
        Popt.lattice.lattice_type, Popt.lattice.lattice_constant, lattice_vec_a, h.ax_show3D_prim, h.ax_show3D_unit);
    
    if     hobj_idx == 1
        FAME_GUI_Info_Set('sphere_radius', num2str(sphere_radius), h)
    elseif hobj_idx == 2
        FAME_GUI_Info_Set('cylinder_radius', num2str(cylinder_radius), h)
    end
end