function Call_Back_rewrite_Lattice_constant(hobj,event,h)
lattice_constant.a=str2num(get(h.edit_constant_a,'string'));
lattice_constant.b=str2num(get(h.edit_constant_b,'string'));
lattice_constant.c=str2num(get(h.edit_constant_c,'string'));
lattice_constant.alpha=str2num(get(h.edit_constant_alpha,'string'));
lattice_constant.beta=str2num(get(h.edit_constant_beta,'string'));
lattice_constant.gamma=str2num(get(h.edit_constant_gamma,'string'));
c = get(h.pop_type,'value');
err=0;
switch c
    case 1
        lattice_type     = 'simple_cubic';
        lattice_constant.b=lattice_constant.a;
        lattice_constant.c=lattice_constant.a;
    case 2
        lattice_type     = 'face_centered_cubic';
        lattice_constant.b=lattice_constant.a;
        lattice_constant.c=lattice_constant.a;
    case 3
        lattice_type     = 'body_centered_cubic';
        lattice_constant.b=lattice_constant.a;
        lattice_constant.c=lattice_constant.a;
    case 4
        lattice_type     = 'hexagonal';
        lattice_constant.b=lattice_constant.a;
    case 5
        lattice_type     = 'rhombohedral';
        lattice_constant.b=lattice_constant.a;
    case 6
        lattice_type     = 'primitive_tetragonal'; 
        lattice_constant.b=lattice_constant.a;
        if lattice_constant.c==lattice_constant.a
            warndlg('  In tetragonal system, lattice constant must satisfying a = b ¡Ú c, £\=£]=£^=£k/2','Warning')
            err=1;
        end
    case 7
        lattice_type     = 'body_centered_tetragonal';
        lattice_constant.b=lattice_constant.a;
        if lattice_constant.c==lattice_constant.a
            warndlg('  In tetragonal system, lattice constant must satisfying a = b ¡Ú c, £\=£]=£^=£k/2','Warning')
            err=1;
        end
    case 8
        lattice_type     = 'primitive_orthorhombic';
        if (lattice_constant.a == lattice_constant.b) || (lattice_constant.a == lattice_constant.c) || (lattice_constant.c == lattice_constant.b)
            warndlg('  In orthorhombic system, lattice constant must satisfying a ¡Ú b ¡Ú c, £\=£]=£^=£k/2','Warning')
            err=1;
        end
    case 9
        lattice_type     = 'face_centered_orthorhombic';
        if (lattice_constant.a == lattice_constant.b) || (lattice_constant.a == lattice_constant.c) || (lattice_constant.c == lattice_constant.b)
            warndlg('  In orthorhombic system, lattice constant must satisfying a ¡Ú b ¡Ú c, £\=£]=£^=£k/2','Warning')
            err=1;
        end
    case 10
        lattice_type     = 'body_centered_orthorhombic';
        if (lattice_constant.a == lattice_constant.b) || (lattice_constant.a == lattice_constant.c) || (lattice_constant.c == lattice_constant.b)
            warndlg('  In orthorhombic system, lattice constant must satisfying a ¡Ú b ¡Ú c, £\=£]=£^=£k/2','Warning')
            err=1;
        end
    case 11
        lattice_type     = 'a_base_centered_orthorhombic';
        if (lattice_constant.a == lattice_constant.b) || (lattice_constant.a == lattice_constant.c) || (lattice_constant.c == lattice_constant.b)
            warndlg('  In orthorhombic system, lattice constant must satisfying a ¡Ú b ¡Ú c, £\=£]=£^=£k/2','Warning')
            err=1;
        end
    case 12
        lattice_type     = 'c_base_centered_orthorhombic';
        if (lattice_constant.a == lattice_constant.b) || (lattice_constant.a == lattice_constant.c) || (lattice_constant.c == lattice_constant.b)
            warndlg('  In orthorhombic system, lattice constant must satisfying a ¡Ú b ¡Ú c, £\=£]=£^=£k/2','Warning')
            err=1;
        end
    case 13
        lattice_type     = 'primitive_monoclinic';
        if (lattice_constant.alpha >= pi/2) || (lattice_constant.a > lattice_constant.c) || (lattice_constant.b > lattice_constant.c)
            warndlg('  In monoclinic system, lattice constant must satisfying a , b ¡Ø c, £\<£k/2, £]=£^=£k/2','Warning')
            err=1;
        end
    case 14
        lattice_type     = 'base_centered_monoclinic';
        if (lattice_constant.alpha >= pi/2) || (lattice_constant.a > lattice_constant.c) || (lattice_constant.b > lattice_constant.c)
            warndlg('  In monoclinic system, lattice constant must satisfying a , b ¡Ø c, £\<£k/2, £]=£^=£k/2','Warning')
            err=1;
        end
    case 15
        lattice_type     = 'triclinic';
end
if err~=1
    set(h.edit_constant_a,'string',num2str(lattice_constant.a));
    set(h.edit_constant_b,'string',num2str(lattice_constant.b));
    set(h.edit_constant_c,'string',num2str(lattice_constant.c));
    set(h.edit_constant_alpha,'string',num2str(lattice_constant.alpha));
    set(h.edit_constant_beta,'string',num2str(lattice_constant.beta));
    set(h.edit_constant_gamma,'string',num2str(lattice_constant.gamma));
    
    grid_num=[16,16,16];
    [ ~, Par_lattice, ~, ~ ] = FAME_Parameter_Generator( grid_num, lattice_type, lattice_constant, []);
    lattice_vec_a    = Par_lattice.lattice_vec_a;
    lattice_constant = Par_lattice.lattice_constant;

    axes(h.ax_display_shape_primitive)
    cla reset
    axes(h.ax_display_shape_unit)
    cla reset
    FAME_Plot_Material_User_Defined(lattice_vec_a(:,1),lattice_vec_a(:,2),lattice_vec_a(:,3), lattice_type, lattice_constant, [], h.ax_display_shape_primitive, 'primitive_cell');
    FAME_Plot_Material_User_Defined(lattice_vec_a(:,1),lattice_vec_a(:,2),lattice_vec_a(:,3), lattice_type, lattice_constant, [], h.ax_display_shape_unit, 'unit_cell');  
    set(  h.clear        , 'callback', {@Call_Back_clear,h}   );
end
end
