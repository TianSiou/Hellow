function Call_Back_element_edit(hobj,event,h)
ele_list = get(h.element_parameter,'string');
if isempty(ele_list) ~= 1
    k        = get(h.element_parameter,'value' );

    set(h.element_edit{1},'string',ele_list(k));
    ele_list = get(h.element_color,'string');
    set(h.element_edit{2},'string',ele_list(k,:));
    ele_list = get(h.element_permitt,'string');
    set(h.element_edit{3},'string',ele_list(k,:));
    ele_list = get(h.element_permeab,'string');
    set(h.element_edit{4},'string',ele_list(k,:));
    ele_list = get(h.element_xi,'string');
    set(h.element_edit{5},'string',ele_list(k,:));
    ele_list = get(h.element_zeta,'string');
    set(h.element_edit{6},'string',ele_list(k,:));
end
end