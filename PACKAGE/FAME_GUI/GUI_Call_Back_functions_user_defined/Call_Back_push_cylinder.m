function Call_Back_push_cylinder(hobj,event,h)
    va   = get(h.element_parameter,'value');
    temp = get(h.list_parameter,'string');
    Script_str;    
    temp = get(h.object_tempinfo_sphere,'string');
    Script_str;

    elt  = get(h.element_edit{1}  ,'string');

    if isempty(elt)~=1
        list= {'[0,0,0]','[1,1,1]','0.04'};
        answer = inputdlg({'Cylinder bot center:','Cylinder top center','Cylinder radius:'},'Cylinder setting',[1,40],list);

        if isempty(answer) ~= 1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %process string from answer( sphere center)
            temp_answer = num2str(eval(['[',answer{1},']']));
            answer{1} = normalize_idx(temp_answer);
            temp_answer = [];
            temp_answer = num2str(eval(['[',answer{2},']']));
            answer{2} = normalize_idx(temp_answer);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            current_material_string = get(h.list_current_material_cylinder,'string');
            number=length(current_material_string)+1;
            current_material_string{end+1} = ['c',num2str(number),'~',char(elt),':[',answer{1},']-[', answer{2},'], R=',answer{3}];
            set(h.list_current_material_cylinder,'string',current_material_string);
            set(h.list_current_material_cylinder,'Value',length(current_material_string));
            %% 
            %grid_nums = str2num(get(h.edit_grid,'string'));
            grid_num = [16,16,16];

            Script_Call_Lattice_Constant
            Script_Lattice_type

%             [ lattice_vec_a, ~, ~, ~, ~ ] = FAME_Parameter_Generator( grid_nums, lattice_type, lattice_constant );
            [ ~, Par_lattice, ~, ~ ] = FAME_Parameter_Generator( grid_num, lattice_type, lattice_constant, []);
            lattice_vec_a    = Par_lattice.lattice_vec_a;
            lattice_constant = Par_lattice.lattice_constant;
        

            eval(['A_temp=[',answer{1},'];']);
            material.parameters{1}.cylinder_bot_centers = A_temp;
            eval(['A_temp=[',answer{2},'];']);
            material.parameters{1}.cylinder_top_centers = A_temp;
            material.parameters{1}.cylinder_radius      = str2double(answer{3});
            temp_c=get(h.element_color,'string');
            eval(['color_temp=[',char(temp_c(va)),'];']);
            material.parameters{1}.color_map            = color_temp/255;

                        %%
            clear temp
            temp = get(h.object_tempinfo_cylinder_material,'string');
            temp{end+1} = ['[',num2str(str2num(temp_c{va})/255),']'];
            set(h.object_tempinfo_cylinder_material,'string',temp);
            clear temp
            temp = get(h.object_tempinfo_cylinder_top,'string');
            temp{end+1} = answer{2};
            set(h.object_tempinfo_cylinder_top,'string',temp);
            clear temp
            temp = get(h.object_tempinfo_cylinder_bot,'string');
            temp{end+1} = answer{1};
            set(h.object_tempinfo_cylinder_bot,'string',temp);
            clear temp
            temp = get(h.object_tempinfo_cylinder_radius,'string');
            temp{end+1} = answer{3};
            set(h.object_tempinfo_cylinder_radius,'string',temp);
            clear temp
            
            FAME_Plot_Material_User_Defined(lattice_vec_a(:,1),lattice_vec_a(:,2),lattice_vec_a(:,3), lattice_type, lattice_constant, material.parameters, h.ax_display_shape_primitive, 'primitive_cell');
            FAME_Plot_Material_User_Defined(lattice_vec_a(:,1),lattice_vec_a(:,2),lattice_vec_a(:,3), lattice_type, lattice_constant, material.parameters, h.ax_display_shape_unit, 'unit_cell');
        end
    end
end