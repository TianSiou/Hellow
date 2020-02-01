function Call_Back_lattice2(hobj,event,h)
    c_shape = get(h.popup_shape,'value');
    switch c_shape
        case {1, 2} %SC SW
            pt_str = 'GXMR';
            path_str = get(h.edit_lattice_show,'string');
            c = fool_check1(pt_str,path_str);
            if c == 0
                set(h.edit_lattice_show,'string','GXMRG');
            end
        case 3 %FCC
            pt_str = 'XULGWK';
            path_str = get(h.edit_lattice_show,'string');
            c = fool_check1(pt_str,path_str);
            if c == 0
                set(h.edit_lattice_show,'string','XULGWKX');
            end
        case 4 %BCC
            pt_str = 'LHPNZQJ';
            path_str = get(h.edit_lattice_show,'string');
            c = fool_check1(pt_str,path_str);
            if c == 0
                set(h.edit_lattice_show,'string','LHPNZQJL');
            end
    end
    
end
function c = fool_check1(pt_str,path_str)
    c = 1;
    for i = 1:length(path_str)
        temp = strfind(pt_str,path_str(i));
        if isempty(temp)
            temp = 0;
        end
        c = c && temp;
    end
    if c == 0
        warndlg('An error input! Please confirm the input path string and uppercase/lowercase.','Warning!');
    end
end