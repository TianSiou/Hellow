function Call_Back_object_delete(hobj,event,h)
graph=0;
num = questdlg('Choose object type','Delete Object','Sphere','Cylinder','Cancel','default');
if isempty(num)~=1
    switch num
        case 'Sphere'
            Sphere_temp = get(h.list_current_material_sphere,'string');
            if isempty(Sphere_temp)~=1
                [sel, ok] = listdlg(...
                        'ListString'    ,Sphere_temp,...
                        'Name'          ,'Choose object to delete',...
                        'OKString'      ,'delete',...
                        'CancelString'  ,'Cancel',...
                        'SelectionMode' ,'Single',...
                        'ListSize'      ,[200 100]);
                if ok==1
                    graph=1;
                    Sphere_material_temp = get(h.object_tempinfo_sphere_material,'string');
                    Sphere_center_temp   = get(h.object_tempinfo_sphere_center,'string');
                    Sphere_radius_temp   = get(h.object_tempinfo_sphere_radius,'string');
                    Sphere_material_temp(sel)=[];
                    Sphere_center_temp(sel)=[];
                    Sphere_radius_temp(sel)=[];
                    Sphere_temp(sel)=[];
                    set(h.object_tempinfo_sphere_material,'string',Sphere_material_temp);
                    set(h.object_tempinfo_sphere_center,'string',Sphere_center_temp);
                    set(h.object_tempinfo_sphere_radius,'string',Sphere_radius_temp);
                    set(h.list_current_material_sphere,'value',length(Sphere_temp));                    
                    set(h.list_current_material_sphere,'string',Sphere_temp);

                end
            end
        case 'Cylinder'
            Cylinder_temp = get(h.list_current_material_cylinder,'string');
            if isempty(Cylinder_temp)~=1
                [sel, ok] = listdlg(...
                        'ListString'    ,Cylinder_temp,...
                        'Name'          ,'Choose object to delete',...
                        'OKString'      ,'delete',...
                        'CancelString'  ,'Cancel',...
                        'SelectionMode' ,'Single',...
                        'ListSize'      ,[200 100]);
                if ok==1
                    graph=1;
                    Cylinder_material_temp = get(h.object_tempinfo_cylinder_material,'string');
                    Cylinder_top_temp   = get(h.object_tempinfo_cylinder_top,'string');
                    Cylinder_bot_temp   = get(h.object_tempinfo_cylinder_bot,'string');
                    Cylinder_radius_temp   = get(h.object_tempinfo_cylinder_radius,'string');
                    Cylinder_material_temp(sel)=[];
                    Cylinder_top_temp(sel)=[];
                    Cylinder_bot_temp(sel)=[];
                    Cylinder_radius_temp(sel)=[];
                    Cylinder_temp(sel)=[];
                    set(h.object_tempinfo_cylinder_material,'string',Cylinder_material_temp);
                    set(h.object_tempinfo_cylinder_top,'string',Cylinder_top_temp);
                    set(h.object_tempinfo_cylinder_bot,'string',Cylinder_bot_temp);
                    set(h.object_tempinfo_cylinder_radius,'string',Cylinder_radius_temp);
                    set(h.list_current_material_cylinder,'value',length(Cylinder_temp)); 
                    set(h.list_current_material_cylinder,'string',Cylinder_temp);
                end
            end
    end
end
%%redraw
if graph==1
    axes(h.ax_display_shape_primitive)
    cla reset
    axes(h.ax_display_shape_unit)
    cla reset
    Script_Call_Lattice_Constant
    Script_Lattice_type
    grid_num = [16,16,16];
    [ ~, Par_lattice, ~, ~ ] = FAME_Parameter_Generator( grid_num, lattice_type, lattice_constant, []);
    lattice_vec_a    = Par_lattice.lattice_vec_a;
    lattice_constant = Par_lattice.lattice_constant;
    temp = get(h.list_parameter,'string');
    Script_str;
    if isempty(get(h.list_current_material_sphere,'string'))~=1
        switch num
            case 'Cylinder'
                Sphere_material_temp = get(h.object_tempinfo_sphere_material,'string');
                Sphere_center_temp   = get(h.object_tempinfo_sphere_center,'string');
                Sphere_radius_temp   = get(h.object_tempinfo_sphere_radius,'string');
        end
        for j = 1:length(Sphere_material_temp)
            material.parameters{j}.color_map = str2num(Sphere_material_temp{j});
            material.parameters{j}.sphere_centers= str2num(Sphere_center_temp{j});
            material.parameters{j}.sphere_radius = str2num(Sphere_radius_temp{j});
            material.parameters{j}.color_map
        end
        FAME_Plot_Material_User_Defined(lattice_vec_a(:,1),lattice_vec_a(:,2),lattice_vec_a(:,3), lattice_type, lattice_constant, material.parameters, h.ax_display_shape_primitive, 'primitive_cell');
        FAME_Plot_Material_User_Defined(lattice_vec_a(:,1),lattice_vec_a(:,2),lattice_vec_a(:,3), lattice_type, lattice_constant, material.parameters, h.ax_display_shape_unit, 'unit_cell');
    end
    if isempty(get(h.list_current_material_cylinder,'string'))~=1
        switch num
            case 'Sphere'
                Cylinder_material_temp = get(h.object_tempinfo_cylinder_material,'string');
                Cylinder_top_temp   = get(h.object_tempinfo_cylinder_top,'string');
                Cylinder_bot_temp   = get(h.object_tempinfo_cylinder_bot,'string');
                Cylinder_radius_temp   = get(h.object_tempinfo_cylinder_radius,'string');
        end
        for j= 1:length(Cylinder_material_temp)
            material.parameters2{j}.color_map = str2num(Cylinder_material_temp{j});
            material.parameters2{j}.cylinder_bot_centers= str2num(Cylinder_bot_temp{j});
            material.parameters2{j}.cylinder_top_centers= str2num(Cylinder_top_temp{j});
            material.parameters2{j}.cylinder_radius = str2num(Cylinder_radius_temp{j});
        end
        FAME_Plot_Material_User_Defined(lattice_vec_a(:,1),lattice_vec_a(:,2),lattice_vec_a(:,3), lattice_type, lattice_constant, material.parameters2, h.ax_display_shape_primitive, 'primitive_cell');
        FAME_Plot_Material_User_Defined(lattice_vec_a(:,1),lattice_vec_a(:,2),lattice_vec_a(:,3), lattice_type, lattice_constant, material.parameters2, h.ax_display_shape_unit, 'unit_cell');
    end
end
end