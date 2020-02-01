function Call_Back_push_sphere(hobj,event,h)
    va        = get(h.element_parameter,'value');
    elt       = get(h.element_edit{1}  ,'string');
    color     = get(h.element_edit{2}  ,'string');
    
if isempty(elt)~=1 && isempty(color)~=1 
    list={'[0,0,0]','0.1'};
    answer = inputdlg({'Sphere center:','Sphere radius:'},'Sphere setting',[1,40],list);
    if isempty(answer) ~= 1 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %process string from answer( sphere center)
        answer{1} = normalize_idx(answer{1});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        current_material_string = get(h.list_current_material_sphere,'string');
        number=length(current_material_string)+1;
        temp=get(h.object_tempinfo_sphere,'string');
        len =length(temp);
        temp{end+1}= ['b',num2str(number),'=[',answer{1},']'];
        current_material_string{end+1} = ['b',num2str(number),':',char(elt),'- [', answer{1},'], R=',answer{2}];
        set(h.list_current_material_sphere,'string',current_material_string);
        set(h.object_tempinfo_sphere,'string',temp);
        set(h.object_tempinfo_sphere,'value',len);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %%temp_save information
        clear temp
        temp = get(h.object_tempinfo_sphere_material,'string');
        temp{end+1} = ['[',num2str(str2num(color{1})/255),']'];
        set(h.object_tempinfo_sphere_material,'string',temp);
        clear temp
        temp = get(h.object_tempinfo_sphere_center,'string');
        temp{end+1} = answer{1};
        set(h.object_tempinfo_sphere_center,'string',temp);
        clear temp
        temp = get(h.object_tempinfo_sphere_radius,'string');
        temp{end+1} = answer{2};
        set(h.object_tempinfo_sphere_radius,'string',temp);
        clear temp
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        grid_num = str2num(get(h.edit_grid,'string'));
        
        Script_Call_Lattice_Constant %Åª¨úlattice constant
        Script_Lattice_type %Åª¨ú lattice type 

%         [ lattice_vec_a, ~, ~, ~, ~ ] = FAME_Parameter_Generator( grid_nums, lattice_type, lattice_constant );
        [ ~, Par_lattice, ~, ~ ] = FAME_Parameter_Generator( grid_num, lattice_type, lattice_constant, []);
        lattice_vec_a    = Par_lattice.lattice_vec_a;
        lattice_constant = Par_lattice.lattice_constant;
        
        temp = get(h.list_parameter,'string');
        Script_str;
        eval(['A_temp=[',answer{1},'];']);
        material.parameters{1}.sphere_centers =A_temp;
        material.parameters{1}.sphere_radius  = str2double(answer{2});

        eval(['color_temp=[',char(color),'];']);
        material.parameters{1}.color_map = color_temp/255;
        FAME_Plot_Material_User_Defined(lattice_vec_a(:,1),lattice_vec_a(:,2),lattice_vec_a(:,3), lattice_type, lattice_constant, material.parameters, h.ax_display_shape_primitive, 'primitive_cell');
        FAME_Plot_Material_User_Defined(lattice_vec_a(:,1),lattice_vec_a(:,2),lattice_vec_a(:,3), lattice_type, lattice_constant, material.parameters, h.ax_display_shape_unit, 'unit_cell');
    end
else
    msgbox('You need to add element information first!', 'Error','error');
end
end