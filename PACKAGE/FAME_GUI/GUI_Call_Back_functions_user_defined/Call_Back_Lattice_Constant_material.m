function Call_Back_Lattice_Constant_material(hobj,event,h)
type_n = get(h.pop_type,'value' ) ;
temp=['material.lattice_constant.a = ',str2double(get(h.edit_constant_a,'string')),';'];
switch type_n
    case 4
        temp=[temp;
            ['material.lattice_constant.c = ',str2double(get(h.edit_constant_c,'string')),';']];
    case 5
        temp=[temp;
            ['material.lattice_constant.c = ',str2double(get(h.edit_constant_c,'string')),';']];
    case 6
        temp=[temp;
            ['material.lattice_constant.c = ',str2double(get(h.edit_constant_c,'string')),';']];
    case 7
        temp=[temp;
            ['material.lattice_constant.c = ',str2double(get(h.edit_constant_c,'string')),';']];
    case 8
        temp=[temp;
            ['material.lattice_constant.b = ',str2double(get(h.edit_constant_b,'string')),';'];
            ['material.lattice_constant.c = ',str2double(get(h.edit_constant_c,'string')),';'] ];
    case 9
        temp=[temp;
            ['material.lattice_constant.b = ',str2double(get(h.edit_constant_b,'string')),';'];
            ['material.lattice_constant.c = ',str2double(get(h.edit_constant_c,'string')),';'] ];
    case 10
        temp=[temp;
            ['material.lattice_constant.b = ',str2double(get(h.edit_constant_b,'string')),';'];
            ['material.lattice_constant.c = ',str2double(get(h.edit_constant_c,'string')),';'] ];
    case 11
        temp=[temp;
            ['material.lattice_constant.b = ',str2double(get(h.edit_constant_b,'string')),';'];
            ['material.lattice_constant.c = ',str2double(get(h.edit_constant_c,'string')),';'] ];
    case 12
        temp=[temp;
            ['material.lattice_constant.b = ',str2double(get(h.edit_constant_b,'string')),';'];
            ['material.lattice_constant.c = ',str2double(get(h.edit_constant_c,'string')),';'] ];
    case 13
        temp=[temp;
            ['material.lattice_constant.b = ',str2double(get(h.edit_constant_b,'string')),';'];
            ['material.lattice_constant.c = ',str2double(get(h.edit_constant_c,'string')),';'];
            ['material.lattice_constant.gamma = ,',str2double(get(h.edit_constant_gamma,'string')),';'] ];
    case 14
        temp=[temp;
            ['material.lattice_constant.b = ',str2double(get(h.edit_constant_b,'string')),';'];
            ['material.lattice_constant.c = ',str2double(get(h.edit_constant_c,'string')),';'];
            ['material.lattice_constant.gamma = ,',str2double(get(h.edit_constant_gamma,'string')),';'] ];
    case 15
        temp=[temp;
            ['material.lattice_constant.b = ',str2double(get(h.edit_constant_b,'string')),';'];
            ['material.lattice_constant.c = ',str2double(get(h.edit_constant_c,'string')),';'];
            ['material.lattice_constant.alpha = ,',str2double(get(h.edit_constant_alpha,'string')),';'];
            ['material.lattice_constant.beta  = ,',str2double(get(h.edit_constant_alpha,'string')),';'];
            ['material.lattice_constant.gamma = ,',str2double(get(h.edit_constant_alpha,'string')),';'] ];
end
set(h.text_Lattice_Constant_material,'string',temp);
