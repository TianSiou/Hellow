function Call_Back_load(hobj,event,h)
    Call_Back_clear(hobj,event,h)
[file_name, path] = uigetfile('*.m','Load file name');
if path~=0
    addpath(path)
    FAME_Folder_Manager
    eval(['U=@(sphere_radius, cylinder_radius) ',sscanf(file_name,'%[^.]'),'(sphere_radius, cylinder_radius);'])
    sphere_radius=.5*ones(100,1); cylinder_radius=.25*ones(100,1);
    U_temp=U(sphere_radius, cylinder_radius);
    sphere_radius=.1*ones(U_temp.material_num,1)*U_temp.lattice_constant.a;
    cylinder_radius=.05*ones(U_temp.material_num,1)*U_temp.lattice_constant.a;
    eval(['material=',sscanf(file_name,'%[^.]'),'(sphere_radius, cylinder_radius);'])
    %% set radius for sphere and cylinder (defalt)
    %%
     n_material = material.material_num;
    %% check lattice type and input lattice constant
    %預設關閉
%     for i = 1:length(h.frame_constant)
%              set(h.frame_constant{i},'visible','off');
%     end
    %讀取lattice type
    set(h.edit_constant_a,'string',num2str(material.lattice_constant.a));
    switch material.lattice_type
                case 'simple_cubic'
                    type = 1;
                case 'face_centered_cubic'
                    type = 2;
                case 'body_centered_cubic'
                    type = 3;
                case 'hexagonal'
                    type = 4;
                    set(h.edit_constant_c,'string',num2str(material.lattice_constant.c));
                case 'rhombohedral'
                    type = 5;
                    set(h.edit_constant_c,'string',num2str(material.lattice_constant.c));
                case 'primitive_tetragonal'
                    type = 6;
                    set(h.edit_constant_c,'string',num2str(material.lattice_constant.c));
                case 'body_centered_tetragonal'
                    type = 7;
                    set(h.edit_constant_c,'string',num2str(material.lattice_constant.c));
                case 'primitive_orthorhombic'
                    type = 8;
                    set(h.edit_constant_b,'string',num2str(material.lattice_constant.b));
                    set(h.edit_constant_c,'string',num2str(material.lattice_constant.c));
                case 'face_centered_orthorhombic'
                    type = 9;
                    set(h.edit_constant_b,'string',num2str(material.lattice_constant.b));
                    set(h.edit_constant_c,'string',num2str(material.lattice_constant.c));
                case 'body_centered_orthorhombic'
                    type = 10;
                    set(h.edit_constant_b,'string',num2str(material.lattice_constant.b));
                    set(h.edit_constant_c,'string',num2str(material.lattice_constant.c));
                case 'a_base_centered_orthorhombic'
                    type = 11;
                    set(h.edit_constant_b,'string',num2str(material.lattice_constant.b));
                    set(h.edit_constant_c,'string',num2str(material.lattice_constant.c));
                case 'c_base_centered_orthorhombic'
                    type = 12;
                    set(h.edit_constant_b,'string',num2str(material.lattice_constant.b));
                    set(h.edit_constant_c,'string',num2str(material.lattice_constant.c));
                case 'primitive_monoclinic'
                    type = 13;
                    set(h.edit_constant_b,'string',num2str(material.lattice_constant.b));
                    set(h.edit_constant_c,'string',num2str(material.lattice_constant.c));
                    set(h.edit_constant_gamma,'string',num2str(material.lattice_constant.gamma));
                case 'a_base_centered_monoclinic'
                    type = 14;
                    set(h.edit_constant_b,'string',num2str(material.lattice_constant.b));
                    set(h.edit_constant_c,'string',num2str(material.lattice_constant.c));
                    set(h.edit_constant_gamma,'string',num2str(material.lattice_constant.gamma));                
                case 'triclinic'
                    type = 15;
                    set(h.edit_constant_b,'string',num2str(material.lattice_constant.b));
                    set(h.edit_constant_c,'string',num2str(material.lattice_constant.c));
                    set(h.edit_constant_alpha,'string',num2str(material.lattice_constant.alpha));
                    set(h.edit_constant_beta,'string',num2str(material.lattice_constant.beta));
                    set(h.edit_constant_gamma,'string',num2str(material.lattice_constant.gamma));
    end
    %選單切換至對應的lattice type，並且打開對應type的lattice constant界面
    set(h.pop_type,'Value',type);
%     set(h.frame_constant{type},'visible','on')
     %% put scalar in to paramter part
     if isfield(material,'scalar') ~=0
      temp = struct2table(material.scalar);
      temp_1 = temp.Properties.VariableNames;
      temp_2 = table2cell(temp);
      n_scalar = length(temp_1) ;
     
   
  
    T_scalar = []; T_scalar_material = [];
    clear temp
    for i = 1 : n_scalar
        temp = [temp_1{i},'=',num2str(temp_2{i}),';'];
        temp_material = ['material.scalar.',temp_1{i},'=',num2str(temp_2{i}),';'];
        T_scalar{end+1} = temp;
        T_scalar_material{end+1} = temp_material;
    end
    
    set(h.list_parameter,'string',T_scalar);
    set(h.list_parameter,'value',n_scalar);
    set(h.list_parameter_material,'string',T_scalar_material);
    set(h.list_parameter_material,'value',n_scalar)
     end
    %% 
    clear temp
    temp_ele = []; temp_color = []; temp_permitt = []; temp_permeab = []; temp_xi = []; temp_zeta = [];
    temp_sphere = []; temp_cylinder =[]; ball_num=0; cylinder_num=0;
    tempinfo_sphere_material = [];      tempinfo_cylinder_material = [];
    tempinfo_sphere_center   = [];      tempinfo_cylinder_bot      = [];       tempinfo_cylinder_top = [];
    tempinfo_sphere_radius   = [];      tempinfo_cylinder_radius     = []; 
    grid_num = str2num(get(h.edit_grid,'string'));
%     [ lattice_vec_a, ~, ~, ~, ~ ] = FAME_Parameter_Generator( grid_nums, material.lattice_type, material.lattice_constant );
    [ ~, Par_lattice, ~, ~ ] = FAME_Parameter_Generator( grid_num, material.lattice_type, material.lattice_constant, []);
    lattice_vec_a    = Par_lattice.lattice_vec_a;
    material.lattice_constant = Par_lattice.lattice_constant;
    for i=1:material.material_num
        temp_ele{end+1} =material.parameters{i}.name;
        temp_color{end+1,1} = ['[',num2str(material.parameters{i}.color_map),']'];
        temp_permitt{end+1} =  13 ;
        temp_permeab{end+1} =   1 ;
        temp_xi{end+1}     = 0 ;
        temp_zeta{end+1}    = 0 ;
        n_sphere = size(material.parameters{i}.sphere_centers,1);
        %%% 檢查球個數、柱子個數
        if isfield(eval(['material.parameters{',num2str(i),'}']),'sphere_centers')==1
            n_sphere = size(material.parameters{i}.sphere_centers,1);
        else
            n_sphere = 0;
        end
        if isfield(eval(['material.parameters{',num2str(i),'}']),'cylinder_top_centers')==1
            n_cylinder = size(material.parameters{i}.cylinder_top_centers,1);
        else
            n_cylinder =0;
        end
        %%%
        material.parameters{i}.color_map=material.parameters{i}.color_map/255;
        for k=1:n_sphere  
            ball_num=ball_num+1;
            temp_str=normalize_idx(num2str(material.parameters{i}.sphere_centers(k,:)));
            temp_sphere{end+1} = ['b',num2str(ball_num),':',...
            material.parameters{i}.name,'-[',temp_str,'], R=',num2str(sphere_radius(i))];
            %%tempinfo
            tempinfo_sphere_material{end+1} = ['[',num2str(material.parameters{i}.color_map),']'];
            tempinfo_sphere_center{end+1}   = ['[',num2str(material.parameters{i}.sphere_centers(k,:)),']'];
            tempinfo_sphere_radius{end+1}   = ['[',num2str(sphere_radius(i)),']'];
            
        end
        for j=1:n_cylinder
            cylinder_num=cylinder_num+1;
            temp_str_bot = normalize_idx(num2str(material.parameters{i}.cylinder_bot_centers(j,:)));
            temp_str_top = normalize_idx(num2str(material.parameters{i}.cylinder_top_centers(j,:)));
            temp_cylinder{end+1} = ['c',num2str(cylinder_num),'~',...
            material.parameters{i}.name,':[',temp_str_bot,...
            ']-[',temp_str_top,'], R=',num2str(cylinder_radius(i))];
            tempinfo_cylinder_material{end+1}     = ['[',num2str(material.parameters{i}.color_map),']'];
            tempinfo_cylinder_bot{end+1} = ['[',num2str(material.parameters{i}.cylinder_bot_centers(j,:)),']'];
            tempinfo_cylinder_top{end+1} = ['[',num2str(material.parameters{i}.cylinder_top_centers(j,:)),']'];
            tempinfo_cylinder_radius{end+1}   = ['[',num2str(cylinder_radius(i)),']'];
            
        end 
    end
    FAME_Plot_Material_User_Defined(lattice_vec_a(:,1),lattice_vec_a(:,2),lattice_vec_a(:,3), material.lattice_type, material.lattice_constant, material.parameters, h.ax_display_shape_unit, 'unit_cell');
    FAME_Plot_Material_User_Defined(lattice_vec_a(:,1),lattice_vec_a(:,2),lattice_vec_a(:,3), material.lattice_type, material.lattice_constant, material.parameters, h.ax_display_shape_primitive, 'primitive_cell');

    set(h.list_current_material_sphere,   'string', temp_sphere);
    set(h.list_current_material_sphere,   'value', size(temp_sphere,1));
    set(h.list_current_material_cylinder,   'string', temp_cylinder);
    set(h.list_current_material_cylinder, 'value', size(temp_cylinder,1));
    set(h.element_parameter , 'string', temp_ele);
    set(h.element_color,          'string', (temp_color));
    set(h.element_permitt,        'string', (temp_permitt));
    set(h.element_permeab,        'string', (temp_permeab));
    set(h.element_xi,             'string', (temp_xi));
    set(h.element_zeta,           'string', (temp_zeta));
    set(h.element_parameter    , 'callback', {@Call_Back_element_edit,h}     );
    set(h.object_tempinfo_sphere_material, 'string' , tempinfo_sphere_material);
    set(h.object_tempinfo_sphere_center, 'string' , tempinfo_sphere_center);
    set(h.object_tempinfo_sphere_radius, 'string' , tempinfo_sphere_radius);
    set(h.object_tempinfo_cylinder_material, 'string' , tempinfo_cylinder_material);
    set(h.object_tempinfo_cylinder_bot, 'string' , tempinfo_cylinder_bot);
    set(h.object_tempinfo_cylinder_top, 'string' , tempinfo_cylinder_top);
    set(h.object_tempinfo_cylinder_radius, 'string' , tempinfo_cylinder_radius);

end

end

