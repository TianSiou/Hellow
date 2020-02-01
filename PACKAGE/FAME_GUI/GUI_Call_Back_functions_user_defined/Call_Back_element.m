function Call_Back_element1(hobj,event,h)
%%使三個選單能選同一個
va = get(h.element_parameter,'value');
% set(h.element_color,'value',va);
% set(h.element_per,'value',va);
ele_list1 = get(h.element_parameter, 'string');
ele_list2 = get(h.element_color,'string');
ele_list3 = get(h.element_per, 'string');
    elt  = char(ele_list3(va));
    eval(['elt=',elt,';']);

S1=['Information of the element ', char(ele_list1(va)) ,':'];
S2=['Color (RGB)  :', '[', char(ele_list2(va)) , ']'  ];
S3=['Permittivity    :', num2str(elt(1)) ];
S4=['Permeability :', num2str(elt(2)) ];

set(h.element_info,'string',{S1;' ';S2;S3;S4});
end