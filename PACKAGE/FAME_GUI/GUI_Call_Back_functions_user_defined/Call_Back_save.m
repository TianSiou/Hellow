function Call_Back_save(hobj,event,h)
%answer = inputdlg({'Sphere center:','Sphere radius:'},'Sphere setting',[1,40],list);
messenge{1}= '%% Basic information';
%% material number
material_list = get(h.element_parameter, 'string');
            n = length(material_list);
messenge{end+1} = ['material.material_num =' ,num2str(n),';'] ;
messenge{end+1} = '';
%% Lattice type
Script_Lattice_type
messenge{end+1} = ['material.lattice_type =''',lattice_type,''';'];
messenge{end+1} = '%% material parameter';
%% Get lattice constant
type_n = get(h.pop_type,'value' ) ;
messenge{end+1}=['material.lattice_constant.a = ',char(get(h.edit_constant_a,'string')),';'];
switch type_n
    case 4
        messenge{end+1}=['material.lattice_constant.c = ',char(get(h.edit_constant_c,'string')),';'];
    case 5
        messenge{end+1}=['material.lattice_constant.c = ',char(get(h.edit_constant_c,'string')),';'];
    case 6
        messenge{end+1}=['material.lattice_constant.c = ',char(get(h.edit_constant_c,'string')),';'];
    case 7
        messenge{end+1}=['material.lattice_constant.c = ',char(get(h.edit_constant_c,'string')),';'];
    case 8
        messenge{end+1}=['material.lattice_constant.b = ',char(get(h.edit_constant_b,'string')),';'];
        messenge{end+1}=['material.lattice_constant.c = ',char(get(h.edit_constant_c,'string')),';'];
    case 9
        messenge{end+1}=['material.lattice_constant.b = ',char(get(h.edit_constant_b,'string')),';'];
        messenge{end+1}=['material.lattice_constant.c = ',char(get(h.edit_constant_c,'string')),';'];
    case 10
        messenge{end+1}=['material.lattice_constant.b = ',char(get(h.edit_constant_b,'string')),';'];
        messenge{end+1}=['material.lattice_constant.c = ',char(get(h.edit_constant_c,'string')),';'];
    case 11
        messenge{end+1}=['material.lattice_constant.b = ',char(get(h.edit_constant_b,'string')),';'];
        messenge{end+1}=['material.lattice_constant.c = ',char(get(h.edit_constant_c,'string')),';'];
    case 12
        messenge{end+1}=['material.lattice_constant.b = ',char(get(h.edit_constant_b,'string')),';'];
        messenge{end+1}=['material.lattice_constant.c = ',char(get(h.edit_constant_c,'string')),';'];
    case 13
        messenge{end+1}=['material.lattice_constant.b = ',char(get(h.edit_constant_b,'string')),';'];
        messenge{end+1}=['material.lattice_constant.c = ',char(get(h.edit_constant_c,'string')),';'];
        messenge{end+1}=['material.lattice_constant.gamma = ,',char(get(h.edit_constant_gamma,'string')),';'];
    case 14
        messenge{end+1}=['material.lattice_constant.b = ',char(get(h.edit_constant_b,'string')),';'];
        messenge{end+1}=['material.lattice_constant.c = ',char(get(h.edit_constant_c,'string')),';'];
        messenge{end+1}=['material.lattice_constant.gamma = ,',char(get(h.edit_constant_gamma,'string')),';'];
    case 15
        messenge{end+1}=['material.lattice_constant.b = ',char(get(h.edit_constant_b,'string')),';'];
        messenge{end+1}=['material.lattice_constant.c = ',char(get(h.edit_constant_c,'string')),';'];
        messenge{end+1}=['material.lattice_constant.alpha = ,',char(get(h.edit_constant_alpha,'string')),';'];
        messenge{end+1}=['material.lattice_constant.beta = ,',char(get(h.edit_constant_beta,'string')),';'];
        messenge{end+1}=['material.lattice_constant.gamma = ,',char(get(h.edit_constant_gamma,'string')),';'];
end
%% parameter list
temp_parameter_normalize = get(h.list_parameter_material, 'string');
for i= 1 : length(temp_parameter_normalize) 
   messenge{end+1} = temp_parameter_normalize{i};
end
%% Foolproof
messenge{end+1} = ['FAME_Material_Locate_Parameter_Foolproof(sphere_radius,cylinder_radius, material.material_num);'];
%% material
Color_map_list = get(h.element_color,'string');
sphere_list = get(h.list_current_material_sphere, 'string');
cylinder_list = get(h.list_current_material_cylinder, 'string'); 
n_sphere = length(sphere_list);
n_cylinder = length(cylinder_list);
temp_parameter = get(h.list_parameter,'string');
% 1. normalize parameter name (sphere_list and cylinder_list)
if isempty(temp_parameter_normalize) ~= 1
    for i=1:max(n_sphere,n_cylinder)
        if i<= n_sphere
            for k=1:length(temp_parameter_normalize) %盢把计临Θ砏絛Α
                temp = char(sphere_list(i));
                temp = strrep(temp,sscanf(char(temp_parameter(k)),'%[^=]') , sscanf(char(temp_parameter_normalize(k)),'%[^=]'));
                sphere_list{i} = temp;                  
            end
        end
        if i<= n_cylinder
            for k=1:length(temp_parameter_normalize) %盢把计临Θ砏絛Α
                temp = char(cylinder_list(i));
                temp = strrep(temp,sscanf(char(temp_parameter(k)),'%[^=]') , sscanf(char(temp_parameter_normalize(k)),'%[^=]'));
                cylinder_list{i} = temp;
            end             
        end
    end
end
% 2. material name
 for i=1:n %材i
     messenge{end+1} = ['%% Shape description for ',char(material_list(i))];
     n_material_sphere = 0;
     n_material_cylinder = 0;
     messenge{end+1} = ['material.parameters{',num2str(i),'}.name =''',char(material_list(i)),''';'];
     messenge{end+1} = ['material.parameters{',num2str(i),'}.color_map=[',char(Color_map_list(i)),'];'];
     temp_sphere_list = [];
     temp_cylinder_bot_list = [];
     temp_cylinder_top_list = [];
     for j=1:max(n_sphere,n_cylinder) %sphere/cylinder listい筁耾材i狥﹁ 
         if j<= n_sphere
              if isempty(intersect( sphere_list{j},material_list{i}))~=1 %耞listい﹚ith
                  temp_2 = strrep(sscanf(sphere_list{j},'%*[^[]%[^]]'),'[','');
                  temp_2 = normalize_idx(temp_2);
                  temp_sphere_list=[temp_sphere_list,temp_2,';'];
                  n_material_sphere=n_material_sphere+1;                
              end     
         end
         if j<= n_cylinder
             if isempty(intersect(cylinder_list{j},material_list{i}))~=1 %耞listい﹚ith
                 % collect bot center
                  temp_2 = sscanf(cylinder_list{j},'%*[^[]%[^-]%*[^[]');
                  temp_2 = normalize_idx(temp_2)
                  temp_cylinder_bot_list=[temp_cylinder_bot_list,temp_2,';'];
                  temp_2 = sscanf(cylinder_list{j},'%*[^-]%*[^[]%[^,]');
                  temp_2 = normalize_idx(temp_2)
                  temp_cylinder_top_list=[temp_cylinder_top_list,temp_2,';'];
                  n_material_cylinder=n_material_cylinder+1;
              end     
         end
     end
     if n_material_sphere~=0
         messenge{end+1} = ['material.parameters{',num2str(i),'}.sphere_centers=[',temp_sphere_list,'];']; 
         messenge{end+1} = ['material.parameters{',num2str(i),'}.sphere_radius=sphere_radius(',num2str(i),')*ones(1,',num2str(n_material_sphere),');']; 
     end
     if n_material_cylinder~=0
         messenge{end+1} = ['material.parameters{',num2str(i),'}.cylinder_bot_centers=[',temp_cylinder_bot_list,'];']; 
         messenge{end+1} = ['material.parameters{',num2str(i),'}.cylinder_top_centers=[',temp_cylinder_top_list,'];'];  
         messenge{end+1} = ['material.parameters{',num2str(i),'}.cylinder_radius=cylinder_radius(',num2str(i),')*ones(1,',num2str(n_material_cylinder),');']; 
     end
 end
[file_name, path] = uiputfile('*.m','Save file name');
if path~=0
    str = [path file_name];
    file = fopen(str,'w');
    title =['function material =',sscanf(file_name,'%[^.]'),'(sphere_radius, cylinder_radius)'];
    fprintf(file,'%s\n',title);
    fprintf(file,'%s\n',messenge{1:end});
    fprintf(file,'%s\n','end');
    fclose(file);
end
end
 

