function Call_Back_element_re_edit(hobj,event,h)
ele_list = get(h.element_parameter,'string');
if isempty(ele_list) ~= 1
    k        = get(h.element_parameter,'value' );

    ele_list = get(h.element_edit{1},'string');
    tp       = get(h.element_parameter,'string' );
    ele_list = strrep(tp,tp(k),char(ele_list));
    set(h.element_parameter,'string',ele_list);
    
    ele_list = get(h.element_edit{2},'string');
    tp       = get(h.element_color,'string' );
    ele_list = strrep(tp,tp(k),char(ele_list));
    set(h.element_color,'string',ele_list);
    
    ele_list = get(h.element_edit{3},'string');
    tp       = get(h.element_permitt,'string' );
    ele_list = strrep(tp,tp(k),char(ele_list));
    set(h.element_permitt,'string',ele_list);
    
    ele_list = get(h.element_edit{4},'string');
    tp       = get(h.element_permeab,'string' );
    ele_list = strrep(tp,tp(k),char(ele_list));
    set(h.element_permeab,'string',ele_list);
    
    ele_list = get(h.element_edit{5},'string');
    tp       = get(h.element_xi,'string' );
    ele_list = strrep(tp,tp(k),char(ele_list));
    set(h.element_xi,'string',ele_list);
    
    ele_list = get(h.element_edit{6},'string');
    tp       = get(h.element_zeta,'string' );
    ele_list = strrep(tp,tp(k),char(ele_list));
    set(h.element_zeta,'string',ele_list);
    
end
end