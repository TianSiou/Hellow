function Call_Back_pop_type(hobj,event,h)
    
    c = get(hobj,'value');
    %% open the corresponding constant penal
%     for i = 1:length(h.frame_constant)
%         if c == i
%             set(h.frame_constant{i},'visible','on');
%         else
%             set(h.frame_constant{i},'visible','off');
%         end
%     end
    %% Plot corresponding primitive cell
    %grid_nums = str2num(get(h.edit_grid,'string'));
    grid_num = [16,16,16];
    switch c
        case 1
            lattice_type     = 'simple_cubic';
            set(h.edit_constant_a,'enable','on');
            set(h.edit_constant_b,'enable','off');
            set(h.edit_constant_c,'enable','off');
            set(h.edit_constant_alpha,'enable','off');
            set(h.edit_constant_beta,'enable','off');
            set(h.edit_constant_gamma,'enable','off');
            lattice_constant.a=1;
            String = ['lattice constant must satisfying a = b = c, £\=£]=£^=£k/2'];
        case 2
            lattice_type     = 'face_centered_cubic';
            set(h.edit_constant_a,'enable','on');
            set(h.edit_constant_b,'enable','off');
            set(h.edit_constant_c,'enable','off');
            set(h.edit_constant_alpha,'enable','off');
            set(h.edit_constant_beta,'enable','off');
            set(h.edit_constant_gamma,'enable','off');
            lattice_constant.a=1;
            String = ['lattice constant must satisfying a = b = c, £\=£]=£^=£k/2'];
        case 3
            lattice_type     = 'body_centered_cubic';
            set(h.edit_constant_a,'enable','on');
            set(h.edit_constant_b,'enable','off');
            set(h.edit_constant_c,'enable','off');
            set(h.edit_constant_alpha,'enable','off');
            set(h.edit_constant_beta,'enable','off');
            set(h.edit_constant_gamma,'enable','off');
            lattice_constant.a=1;
            String = ['lattice constant must satisfying a = b = c, £\=£]=£^=£k/2'];
        case 4
            lattice_type     = 'hexagonal';
            set(h.edit_constant_a,'enable','on');
            set(h.edit_constant_b,'enable','off');
            set(h.edit_constant_c,'enable','on');
            set(h.edit_constant_alpha,'enable','off');
            set(h.edit_constant_beta,'enable','off');
            set(h.edit_constant_gamma,'enable','off');
            lattice_constant.a=1;
            lattice_constant.c=1;
            String = ['lattice constant must satisfying a = b, £\=£]=£k/2, £^=2£k/3'];
        case 5
            lattice_type     = 'rhombohedral';
            set(h.edit_constant_a,'enable','on');
            set(h.edit_constant_b,'enable','off');
            set(h.edit_constant_c,'enable','on');
            set(h.edit_constant_alpha,'enable','off');
            set(h.edit_constant_beta,'enable','off');
            set(h.edit_constant_gamma,'enable','off');
            lattice_constant.a=1;
            lattice_constant.c=1;
            String = ['lattice constant must satisfying a = b, £\=£]=£k/2, £^=2£k/3'];
        case 6
            lattice_type     = 'primitive_tetragonal';
            set(h.edit_constant_a,'enable','on');
            set(h.edit_constant_b,'enable','off');
            set(h.edit_constant_c,'enable','on');
            set(h.edit_constant_alpha,'enable','off');
            set(h.edit_constant_beta,'enable','off');
            set(h.edit_constant_gamma,'enable','off');
            lattice_constant.a=1;
            lattice_constant.c=2;   
            String = ['lattice constant must satisfying a = b ¡Ú c, £\=£]=£^=£k/2'];

        case 7
            lattice_type     = 'body_centered_tetragonal';
            set(h.edit_constant_a,'enable','on');
            set(h.edit_constant_b,'enable','off');
            set(h.edit_constant_c,'enable','on');
            set(h.edit_constant_alpha,'enable','off');
            set(h.edit_constant_beta,'enable','off');
            set(h.edit_constant_gamma,'enable','off');
            lattice_constant.a=1;
            lattice_constant.c=2;        
            String = ['lattice constant must satisfying a = b ¡Ú c, £\=£]=£^=£k/2'];
        case 8
            lattice_type     = 'primitive_orthorhombic';
            set(h.edit_constant_a,'enable','on');
            set(h.edit_constant_b,'enable','on');
            set(h.edit_constant_c,'enable','on');
            set(h.edit_constant_alpha,'enable','off');
            set(h.edit_constant_beta,'enable','off');
            set(h.edit_constant_gamma,'enable','off');
            lattice_constant.a=1;
            lattice_constant.b=2;
            lattice_constant.c=3;  
            String = ['lattice constant must satisfying a ¡Ú b ¡Ú c, £\=£]=£^=£k/2'];
        case 9
            lattice_type     = 'face_centered_orthorhombic';
            set(h.edit_constant_a,'enable','on');
            set(h.edit_constant_b,'enable','on');
            set(h.edit_constant_c,'enable','on');
            set(h.edit_constant_alpha,'enable','off');
            set(h.edit_constant_beta,'enable','off');
            set(h.edit_constant_gamma,'enable','off');
            lattice_constant.a=1;
            lattice_constant.b=2;
            lattice_constant.c=3;
            String = ['lattice constant must satisfying a ¡Ú b ¡Ú c, £\=£]=£^=£k/2'];
        case 10
            lattice_type     = 'body_centered_orthorhombic';
            set(h.edit_constant_a,'enable','on');
            set(h.edit_constant_b,'enable','on');
            set(h.edit_constant_c,'enable','on');
            set(h.edit_constant_alpha,'enable','off');
            set(h.edit_constant_beta,'enable','off');
            set(h.edit_constant_gamma,'enable','off');
            lattice_constant.a=1;
            lattice_constant.b=2;
            lattice_constant.c=3;
            String = ['lattice constant must satisfying a ¡Ú b ¡Ú c, £\=£]=£^=£k/2'];
        case 11
            lattice_type     = 'a_base_centered_orthorhombic';
            set(h.edit_constant_a,'enable','on');
            set(h.edit_constant_b,'enable','on');
            set(h.edit_constant_c,'enable','on');
            set(h.edit_constant_alpha,'enable','off');
            set(h.edit_constant_beta,'enable','off');
            set(h.edit_constant_gamma,'enable','off');
            lattice_constant.a=1;
            lattice_constant.b=2;
            lattice_constant.c=3; 
            String = ['lattice constant must satisfying a ¡Ú b ¡Ú c, £\=£]=£^=£k/2'];
        case 12
            lattice_type     = 'c_base_centered_orthorhombic';
            set(h.edit_constant_a,'enable','on');
            set(h.edit_constant_b,'enable','on');
            set(h.edit_constant_c,'enable','on');
            set(h.edit_constant_alpha,'enable','off');
            set(h.edit_constant_beta,'enable','off');
            set(h.edit_constant_gamma,'enable','off');
            lattice_constant.a=1;
            lattice_constant.b=2;
            lattice_constant.c=3;
            String = ['lattice constant must satisfying a ¡Ú b ¡Ú c, £\=£]=£^=£k/2'];
        case 13
            lattice_type     = 'primitive_monoclinic';
            set(h.edit_constant_a,'enable','on');
            set(h.edit_constant_b,'enable','on');
            set(h.edit_constant_c,'enable','on');
            set(h.edit_constant_alpha,'enable','on');
            set(h.edit_constant_beta,'enable','off');
            set(h.edit_constant_gamma,'enable','off');
            lattice_constant.a=1;
            lattice_constant.b=2;
            lattice_constant.c=4;
            lattice_constant.alpha=pi/3;
            String = ['lattice constant must satisfying a, b¡Øc, £\<£k/2, £]=£^=£k/2'];
        case 14
            lattice_type     = 'base_centered_monoclinic';
            set(h.edit_constant_a,'enable','on');
            set(h.edit_constant_b,'enable','on');
            set(h.edit_constant_c,'enable','on');
            set(h.edit_constant_alpha,'enable','on');
            set(h.edit_constant_beta,'enable','off');
            set(h.edit_constant_gamma,'enable','off'); 
            lattice_constant.a=1;
            lattice_constant.b=2;
            lattice_constant.c=4;
            lattice_constant.alpha=pi/3;  
            String = ['lattice constant must satisfying a, b¡Øc, £\<£k/2, £]=£^=£k/2'];
        case 15
            lattice_type     = 'triclinic';
            set(h.edit_constant_a,'enable','on');
            set(h.edit_constant_b,'enable','on');
            set(h.edit_constant_c,'enable','on');
            set(h.edit_constant_alpha,'enable','on');
            set(h.edit_constant_beta,'enable','on');
            set(h.edit_constant_gamma,'enable','on'); 
            lattice_constant.a=1;
            lattice_constant.b=2;
            lattice_constant.c=3;
            lattice_constant.alpha=pi/2;
            lattice_constant.beta=pi/3;
            lattice_constant.gamma=pi/4;   
            String = [];
    end
    set(h.text_information,'string',String);
    [ ~, Par_lattice, ~, ~ ] = FAME_Parameter_Generator( grid_num, lattice_type, lattice_constant, []);
    lattice_vec_a    = Par_lattice.lattice_vec_a;
    lattice_constant = Par_lattice.lattice_constant;
    set(h.edit_constant_a,'string',num2str(lattice_constant.a));
    set(h.edit_constant_b,'string',num2str(lattice_constant.b));
    set(h.edit_constant_c,'string',num2str(lattice_constant.c));
    set(h.edit_constant_alpha,'string',num2str(lattice_constant.alpha));
    set(h.edit_constant_beta,'string',num2str(lattice_constant.beta));
    set(h.edit_constant_gamma,'string',num2str(lattice_constant.gamma));

    axes(h.ax_display_shape_primitive)
    cla reset
    axes(h.ax_display_shape_unit)
    cla reset
    FAME_Plot_Material_User_Defined(lattice_vec_a(:,1),lattice_vec_a(:,2),lattice_vec_a(:,3), lattice_type, lattice_constant, [], h.ax_display_shape_primitive, 'primitive_cell');
    FAME_Plot_Material_User_Defined(lattice_vec_a(:,1),lattice_vec_a(:,2),lattice_vec_a(:,3), lattice_type, lattice_constant, [], h.ax_display_shape_unit, 'unit_cell');  
    
    %% Clear current shape list
    set(  h.clear        , 'callback', {@Call_Back_clear,h}   );
end