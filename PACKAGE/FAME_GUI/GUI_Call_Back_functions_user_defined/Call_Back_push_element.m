function Call_Back_push_element(hobj,event,h)
    list = {'A','[0,0,0]','13','1','0','0'};
    answer = inputdlg({'element name:','color(0~255):','permittivity','permeability','Electromagnetic parameter ξ :','Electromagnetic parameter ζ :'},'Add parameter',[1,40],list);
%    if length(answer{1})>0 && length(answer{2})>0
%    if isempty(answer{1})~=1 && isempty(answer{2})~=1 && isempty(answer{3})~=1 && isempty(answer{4})~=1
     if isempty(answer)~=1
         if isempty(answer{1})~=1 && isempty(answer{2})~=1 && isempty(answer{3})~=1 && isempty(answer{4})~=1
        ele_list = get(h.element_parameter,'string');
        ele_list{end+1} = [answer{1}];
        set(h.element_parameter,'string',ele_list);
        set(h.element_parameter,'value',length(ele_list));
        
        %% 填入元素資訊(暫存)
        ele_list = get(h.element_color,'string');
        ele_list{end+1} = [answer{2}];
        set(h.element_color  ,'string',ele_list);  %color 
        ele_list = get(h.element_permitt,'string');
        ele_list{end+1} = [answer{3}];
        set(h.element_permitt,'string',ele_list);
        ele_list = get(h.element_permeab,'string');
        ele_list{end+1} = [answer{4}];
        set(h.element_permeab,'string',ele_list);
        ele_list = get(h.element_xi,'string');
        ele_list{end+1} = [answer{5}];
        set(h.element_xi     ,'string',ele_list);
        ele_list = get(h.element_zeta,'string');
        ele_list{end+1} = [answer{6}];
        set(h.element_zeta   ,'string',ele_list);
        
        Call_Back_element_edit(hobj,event,h)
         else
           msgbox({'It is not enough information for element!';'Please re-key it.'}, 'Warn','warn');
         end
    end
%         temp = get(h.element_parameter,'string');
%         Script_str
end